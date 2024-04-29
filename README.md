
# pltree
基于持久化内存的学习索引

pltree-local,保留了内部节点的局部搜索

pltree-n,关闭了块元数据的加速

pltree-ff,将溢出缓存替换为fast&fair


# 编译
mkdir build
cd build
cmake ..
make -j 14
# 运行
运行时执行命令：sudo sh test.sh
test.sh中参数说明：
./run.sh /home/zzg/data/uniform_dense_200M_uint64 binary uint64_t 100000000 200000000 mixed 0 1 0 1 >>Pal.txt

/home/zzg/data/uniform_dense_200M_uint64：keys_file 数据集路径

binary：keys_file_type

uint64_t：keys_type

100000000：workload_keys

200000000：total_num_keys

mixed ：operation （mixed ，range，erase，insert，search，update）

第10个参数0:设置写入比例：insert_frac 插入操作:1，搜索操作:0，读写混合(0-1)

run.sh中参数说明：
设置索引初始化加载数据量：  --init_num_keys=100000000 \
倾斜负载设置：
--lookup_distribution=zipf \
--theta=0.99 \

