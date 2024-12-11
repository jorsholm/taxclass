# download, unzip, and build protaxA
git clone https://github.com/psomervuo/protaxA.git
unzip protaxA/c.zip
unzip protaxA/c2.zip
unzip protaxA/scripts.zip
cd c
make CFLAGS="-march=native -mtune=native -O3"
cp create_xdata_best2 trainclassify_best2 classify_best2 ../scripts
cd ../c2
rm *.o classify_info classify_rseq classify_v1 classify_v2
make CFLAGS="-march=native -mtune=native -O3"
