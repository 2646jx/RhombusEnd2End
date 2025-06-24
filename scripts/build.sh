. scripts/common.sh

for deps in eigen3 emp-ot emp-tool hexl SEAL-3.7
do
  if [ ! -d $BUILD_DIR/include/$deps ] 
  then
 echo -e "${RED}$deps${NC} seems absent in ${BUILD_DIR}/include/, please re-run scripts/build-deps.sh"
 exit 1
  fi
done

for deps in zstd.h 
do
  if [ ! -f $BUILD_DIR/include/$deps ] 
  then
 echo -e "${RED}$deps${NC} seems absent in ${BUILD_DIR}/include/, please re-run scripts/build-deps.sh"
 exit 1
  fi
done

cd $BUILD_DIR/
cmake .. -DCMAKE_BUILD_TYPE=Release -DSCI_BUILD_NETWORKS=ON -DSCI_BUILD_TESTS=OFF -DOPENSSL_ROOT_DIR=/usr/local -DCMAKE_PREFIX_PATH=$BUILD_DIR -DUSE_APPROX_RESHARE=ON

make resnet50-cheetah -j4 
make resnet50-rhombus -j4

# rhombus module test
make rhombus_matmul
make rhombus_matvec