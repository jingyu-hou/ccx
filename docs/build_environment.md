# Build Environment

## Current Status

WSL Ubuntu has not been installed yet on this machine. Enabling WSL requires an elevated Windows PowerShell session.

Observed host system:

- Windows: Windows 10 Home China 2009
- WSL: not enabled / not installed
- Required action: run the commands below as Administrator, then reboot if prompted.

```powershell
wsl --install -d Ubuntu
```

If `wsl --install` still reports that WSL is unavailable, enable the required Windows features first:

```powershell
dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart
wsl --set-default-version 2
wsl --install -d Ubuntu
```

## Target Build Environment

- Operating system: WSL Ubuntu xx, pending installation
- gcc: xx, install with `sudo apt install build-essential`
- gfortran: xx, install with `sudo apt install gfortran`
- make: xx, install with `sudo apt install make`
- cmake: xx, install with `sudo apt install cmake`
- BLAS/LAPACK: system apt package versions, install with `sudo apt install libblas-dev liblapack-dev`
- ARPACK: system apt package or source build, prefer `sudo apt install libarpack2-dev` unless project compatibility requires source build
- SPOOLES: source build or existing Linux static library

## Version Commands After Installation

Run these inside Ubuntu:

```bash
lsb_release -a
gcc --version
gfortran --version
make --version
cmake --version
apt-cache policy libblas-dev liblapack-dev libarpack2-dev
```

## Project Notes

The current project Makefiles are GNU/Linux-oriented and expect tools such as `gcc`, `gfortran`, `make`, `ar`, `cp`, `rm`, `pwd`, and `perl`.

The existing ARPACK/LAPACK static libraries under `D:\ccx` are Linux ELF archives, so WSL/Linux is the correct build target. Cygwin is installed and has GCC/GFortran, but Cygwin object files are COFF format and are not ABI-compatible with the existing Linux ELF static libraries.
