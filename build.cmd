@Echo Off
if not exist build mkdir build
cd build
if exist msmpi.dll del /s /q msmpi.dll
cmake .. -G "MinGW Makefiles"
cmake --build .
cd ..
pause