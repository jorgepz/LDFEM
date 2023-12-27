
cd src

echo "Creting objects for ldfem:"

# compilo
gfortran -mcmodel=large -c -g depsmodule.f90
echo "  depsmodule."
gfortran -mcmodel=large -c -g usefull.f90
echo "  usefull."
gfortran -mcmodel=large -c -g declarmodule.f90
echo "  declarmodule."
gfortran -mcmodel=large -c -g materialmodule.f90 declarmodule.f90 depsmodule.f90
echo "  materialmodule."
gfortran -mcmodel=large -c -g vtkmodule.f90 declarmodule.f90
echo "  vtkmodule."
gfortran -mcmodel=large -c -g elemmodule.f90 declarmodule.f90 usefull.f90 depsmodule.f90
echo "  elemmodule."
gfortran -mcmodel=large -c -g boundarymodule.f90 declarmodule.f90 usefull.f90 depsmodule.f90
echo "  boundarymodule."
gfortran -mcmodel=large -c -g ldfem.f90 vtkmodule.f90 elemmodule.f90 declarmodule.f90 usefull.f90 \
  boundarymodule.f90 materialmodule.f90 depsmodule.f90
echo "  ldfem."
echo "OK."

echo "Compilling ldfem.lnx ."

gfortran -mcmodel=large -o ldfem.lnx ldfem.o vtkmodule.o elemmodule.o declarmodule.o usefull.o \
  boundarymodule.o  materialmodule.o depsmodule.o \
  -llapack -lblas
  #~ -L$HOME/libf77/$ARCH -llapack -lblas

echo "OK."

echo "cleaning..."
rm *.o *.mod

mv ldfem.lnx ../

cd ..

echo "objects deleted and ldfem moved to ../."
