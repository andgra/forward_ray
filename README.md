# Сборка проекта

1. Установить необходимый софт:
   - MinGw [https://sourceforge.net/projects/mingw-w64/](https://sourceforge.net/projects/mingw-w64/). При установке оставлять параметры по умолчанию. После установки путь ```C:\Program Files (x86)\mingw-w64\i686-8.1.0-posix-dwarf-rt_v6-rev0\mingw32\bin``` добавить в переменную ```PATH``` системы.
   - MPI [https://www.microsoft.com/en-us/download/details.aspx?id=57467](https://www.microsoft.com/en-us/download/details.aspx?id=57467). Для сборки устанавливать ```SDK``` **и** ```setup```. Для использования устанавливать ```setup``` на конечном компьютере (если оставлять переменную ```INCLUDE_MPI_DLL``` в ```false``` в ```CMakeLists.txt```).
   - Cmake [https://cmake.org/download/](https://cmake.org/download/). Выбирать ```Windows win64-x64 Installer```.
2. Если билдить под компьютер без MPI (development mode), то в файле ```CMakeLists.txt``` переменную ```INCLUDE_MPI_DLL``` установить в ```true```, иначе - оставлять ```false``` (может быть конфликт версий, если всегда включать dll в билд). Этот ```dll``` копируется в билд из ```src/lib/msmpi.dll```
3. Запустить ```build.cmd```

## Если есть CLion (ide для c/c++)

3 пункт (запуск сборки) заменяется встроенным инструментом сборки из ide.

# Запуск программы

- В файле ```run.cmd``` или ```mpi.cmd``` поменять последний параметр (где указан путь к моделям в папку ```data```) на нужную модель.
- Для использования технологии mpi запускать ```mpi.cmd``` (с нужным числом потоков - параметр ```-n``` в редактировании файла)
- Для запуска в единой памяти использовать ```run.cmd```
