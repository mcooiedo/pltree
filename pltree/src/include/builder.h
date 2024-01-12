
#ifndef PLTREE_BUILDER_H
#define PLTREE_BUILDER_H

#include <cassert>
#include <limits>
#include <cmath>
#include <vector>
#include "common.h"
#include "allocator.h"

namespace plt
{

    //const double scale=1.3;
    template <class KeyType>//,class Alloc=my_alloc::PMAllocator
    class Builder
    {
    public:
        Builder(KeyType min_key, KeyType max_key, size_t max_error/*,Alloc &alloc*/)
            : min_key_(min_key),
              max_key_(max_key),
              max_error_(max_error),
              curr_num_keys_(0),
              curr_num_distinct_keys_(0),
              prev_key_(min_key),
              prev_position_(0),
              start_(true)//,alloc(alloc)

        {
        }

        // Adds a key. Assumes that keys are stored in a dense array.
        void AddKey(KeyType key)
        {
            if (curr_num_keys_ == 0)
            {
                SetNextSplineStartPoint(key, 0);
                AddKey(key, /*position=*/0);
                return;
            }
            AddKey(key, (prev_position_ + 1));
        }

        // Finalizes the construction.
        void Finalize()
        {
            // Last key needs to be equal to `max_key_`.
            assert(curr_num_keys_ == 0 || prev_key_ == max_key_);

            // Ensure that `prev_key_` (== `max_key_`) is last key on spline.
            if (curr_num_keys_ > 0 && (segments_message_.empty() || segments_message_.back().key != prev_key_))
            {
                AddSegmentMessage();
                //                AddNotFullSegmentMessage();
            }
        }

        inline  std::vector<SegmentMessage<KeyType>> &get_segments_message()
        {
            return segments_message_;
        }

    private:
        void AddKey(KeyType key, size_t position)
        {
            assert(key >= min_key_ && key <= max_key_);
            // Keys need to be monotonically increasing.
            assert(key >= prev_key_);
            // Positions need to be strictly monotonically increasing.
            assert(position == 0 || position > prev_position_);

            PossiblyAddKeyToSpline(key, position);

            ++curr_num_keys_;
            prev_key_ = key;
            prev_position_ = position;
        }

        void AddSegmentMessage()
        {
             size_t start = spline_start_point_.y, len = curr_num_keys_;
             float slope = (prev_point_.x == spline_start_point_.x) ? 0 : (prev_point_.y - spline_start_point_.y) / (prev_point_.x - spline_start_point_.x);
            segments_message_.push_back({prev_point_.x, start, len, slope});
        }

        //        void AddFullSegmentMessage() {
        //            size_t start = spline_start_point_.y, len = curr_num_keys_;
        //            double slope = (prev_point_.x == spline_start_point_.x) ? 0 : (prev_point_.y - spline_start_point_.y) / (prev_point_.x - spline_start_point_.x);
        //            segments_message_.push_back({prev_point_.x, start, len, slope, true});
        //        }
        //
        //        void AddNotFullSegmentMessage() {
        //            size_t start = spline_start_point_.y, len = curr_num_keys_;
        //            double slope = (prev_point_.x == spline_start_point_.x) ? 0 : (prev_point_.y - spline_start_point_.y) / (prev_point_.x - spline_start_point_.x);
        //            segments_message_.push_back({prev_point_.x, start, len, slope, false});
        //        }

        enum Orientation
        {
            Collinear,
            CW,
            CCW
        };
        static constexpr double precision = std::numeric_limits<double>::epsilon();

        static Orientation ComputeOrientation(const double dx1, const double dy1, const double dx2, const double dy2)
        {
            const double expr = std::fma(dy1, dx2, -std::fma(dy2, dx1, 0));
            if (expr > precision)
                return Orientation::CW;
            else if (expr < -precision)
                return Orientation::CCW;
            return Orientation::Collinear;
        };

        void SetUpperLimit(KeyType key, double position) { upper_limit_ = {key, position}; }
        void SetLowerLimit(KeyType key, double position) { lower_limit_ = {key, position}; }
        void SetPreviousCdfPoint(KeyType key, double position) { prev_point_ = {key, position}; }

        void SetNextSplineStartPoint(KeyType key, double position)
        {
            spline_start_point_ = {key, position};
        }

        // Implementation is based on `GreedySplineCorridor` from:
        // T. Neumann and S. Michel. Smooth interpolating histograms with error guarantees. [BNCOD'08]
        void PossiblyAddKeyToSpline(KeyType key, double position)
        {
            if (start_)
            {
                // Add first CDF point to spline.
                ++curr_num_distinct_keys_;
                SetPreviousCdfPoint(key, position);
                start_ = false;
                return;
            }

            if (key == prev_key_)
            {
                // No new CDF point if the key didn't change.
                return;
            }

            // New CDF point.
            ++curr_num_distinct_keys_;

            if (curr_num_distinct_keys_ == 2)
            {
                // Initialize `upper_limit_` and `lower_limit_` using the second CDF point.
                SetUpperLimit(key, position + max_error_);
                SetLowerLimit(key, (position < max_error_) ? 0 : position - max_error_);
                SetPreviousCdfPoint(key, position);
                return;
            }

            // `B` in algorithm.
            const Coord<KeyType> &last = spline_start_point_;

            // Compute current `upper_y` and `lower_y`.
            const double upper_y = position + max_error_;
            const double lower_y = (position < max_error_) ? 0 : position - max_error_;

            // Compute differences.
            assert(upper_limit_.x >= last.x);
            assert(lower_limit_.x >= last.x);
            assert(key >= last.x);
            const double upper_limit_x_diff = upper_limit_.x - last.x;
            const double lower_limit_x_diff = lower_limit_.x - last.x;
            const double x_diff = key - last.x;

            assert(upper_limit_.y >= last.y);
            assert(position >= last.y);
            const double upper_limit_y_diff = upper_limit_.y - last.y;
            const double lower_limit_y_diff = lower_limit_.y - last.y;
            const double y_diff = position - last.y;

            // `prev_point_` is the previous point on the CDF and the next candidate to be added to the spline.
            // Hence, it should be different from the `last` point on the spline.
            assert(prev_point_.x != last.x);

            // Do we cut the error corridor?
            if ((ComputeOrientation(upper_limit_x_diff, upper_limit_y_diff, x_diff, y_diff) != Orientation::CW) || (ComputeOrientation(lower_limit_x_diff, lower_limit_y_diff, x_diff, y_diff) != Orientation::CCW))
            {
                // Add previous CDF point to spline.
                //                AddFullSegmentMessage();
                AddSegmentMessage();

                curr_num_keys_ = 0; // Add 1 outside the function
                curr_num_distinct_keys_ = 1;
                SetNextSplineStartPoint(key, position);
            }
            else
            {
                assert(upper_y >= last.y);
                const double upper_y_diff = upper_y - last.y;
                if (ComputeOrientation(upper_limit_x_diff, upper_limit_y_diff, x_diff, upper_y_diff) == Orientation::CW)
                {
                    SetUpperLimit(key, upper_y);
                }

                const double lower_y_diff = lower_y - last.y;
                if (ComputeOrientation(lower_limit_x_diff, lower_limit_y_diff, x_diff, lower_y_diff) == Orientation::CCW)
                {
                    SetLowerLimit(key, lower_y);
                }
            }

            SetPreviousCdfPoint(key, position);
        }

        const KeyType min_key_;
        const KeyType max_key_;
        const size_t max_error_;
        //Alloc &alloc;

        std::vector<SegmentMessage<KeyType>> segments_message_;

        size_t curr_num_keys_;
        size_t curr_num_distinct_keys_;
        KeyType prev_key_;
        size_t prev_position_;

        // Current upper and lower limits on the error corridor of the spline.
        Coord<KeyType> upper_limit_;
        Coord<KeyType> lower_limit_;

        // Previous CDF point.
        Coord<KeyType> prev_point_;
        Coord<KeyType> spline_start_point_;

        bool start_;
    };

} // 

#endif 