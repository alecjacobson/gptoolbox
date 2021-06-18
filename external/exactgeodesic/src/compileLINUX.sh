g++ -shared -fPIC -I/opt/local/include/ -I/home/bqrosen/boost_1_63_0 -I./ geodesic_matlab_api.cpp -o geodesic.so ;
cp geodesic.so ../geodesic_release.so
