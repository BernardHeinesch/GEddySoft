@ECHO OFF

pushd %~dp0

REM Command file for Sphinx documentation

if "%SPHINXBUILD%" == "" (
	set SPHINXBUILD=sphinx-build
)
set SOURCEDIR=source
set BUILDDIR=build

where %SPHINXBUILD% >NUL 2>NUL
if errorlevel 1 (
	echo.
	echo Error: '%SPHINXBUILD%' not found. Please install Sphinx, for example: pip install -r requirements.txt, and ensure it is on PATH.
	popd
	exit /b 1
)

%SPHINXBUILD% -M html %SOURCEDIR% %BUILDDIR% %SPHINXOPTS% %O%
if errorlevel 1 (
	popd
	exit /b 1
)

popd

echo.
echo Documentation built in %BUILDDIR%\html
echo Open %BUILDDIR%\html\index.html in your browser to view it
