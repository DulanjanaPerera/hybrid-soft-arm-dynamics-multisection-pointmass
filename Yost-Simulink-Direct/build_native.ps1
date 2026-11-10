param(
    [string]$MatlabRoot = 'C:\Program Files\MATLAB\R2025a',
    [string]$Compiler = 'C:\mingw64\bin\g++.exe'
)
# Optional standalone build when MATLAB batch startup is unavailable.
# Uses MATLAB's installed headers/import libraries; never starts or connects to devices.
$ErrorActionPreference = 'Stop'
$project = $PSScriptRoot
$common = @('-std=c++17','-O2','-fwrapv','-fexceptions','-fno-omit-frame-pointer')
& $Compiler @common '-static' (Join-Path $project 'test_parser.cpp') '-o' (Join-Path $project 'test_parser.exe')
if ($LASTEXITCODE -ne 0) { throw 'Parser test compilation failed.' }
& (Join-Path $project 'test_parser.exe')
if ($LASTEXITCODE -ne 0) { throw 'Parser tests failed.' }
$mexArgs = $common + @('-shared','-static','-DMATLAB_MEX_FILE','-DMX_COMPAT_64',
    "-I$MatlabRoot\extern\include", "-I$MatlabRoot\simulink\include",
    (Join-Path $project 'yost_imu_sfun.cpp'), "$MatlabRoot\extern\version\cpp_mexapi_version.cpp",
    '-Wl,--no-undefined', "-Wl,$MatlabRoot\extern\lib\win64\mingw64\mexFunction.def",
    "-Wl,$MatlabRoot\extern\lib\win64\mingw64\exportsmexfileversion.def",
    "-L$MatlabRoot\extern\lib\win64\mingw64",'-llibmex','-llibmx','-llibmat',
    '-o',(Join-Path $project 'yost_imu_sfun.mexw64'))
& $Compiler @mexArgs
if ($LASTEXITCODE -ne 0) { throw 'MEX compilation failed.' }
Write-Output 'Compiled MEX successfully. MATLAB loading and Simulink execution still require validate_yost.'
