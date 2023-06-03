@echo off
setlocal enabledelayedexpansion

set "search=Inverse"
set "replace=TestPoint1_"

for %%F in (*.*) do (
    set "filename=%%~nxF"
    set "newfilename=!filename:%search%=%replace%!"
    ren "%%F" "!newfilename!"
)

echo 文件名替换完成。
pause
