#Author : Chayan Kumar Saha
#Contact : chayan.kumar@umu.se

#!/bin/sh

if [ "$(uname)" = "Darwin" ]; then
	if [ -f ~/.bash_profile ]; then
		if ! [ -x "$(command -v wget)" ]; then
			/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
			brew install wget
			exit 1
		fi
		if [ -x "$(command -v wget)" ]; then
			wget https://repo.anaconda.com/archive/Anaconda3-5.2.0-MacOSX-x86_64.sh
			bash ./Anaconda3-5.2.0-MacOSX-x86_64.sh
			source ~/.bash_profile
		fi
	fi
	if ! [ -f ~/.bash_profile ]; then
		touch ~/.bash_profile
		if ! [ -x "$(command -v wget)" ]; then
			/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
			brew install wget
			exit 1
		fi
		if [ -x "$(command -v wget)" ]; then
			wget https://repo.anaconda.com/archive/Anaconda3-5.2.0-MacOSX-x86_64.sh
			bash ./Anaconda3-5.2.0-MacOSX-x86_64.sh
			source ~/.bash_profile
		fi
	fi
fi

if [ "$(uname)" = "Linux" ]; then
	if [ "$(uname -m)" = "x86_64" ]; then
		if [ -f ~/.bashrc ]; then
			if ! [ -x "$(command -v wget)" ]; then
				echo 'Kindly install "Wget" and run this script again, thanks!!' >&2
				exit 1
			fi
			if [ -x "$(command -v wget)" ]; then
				wget https://repo.anaconda.com/archive/Anaconda3-5.2.0-Linux-x86_64.sh
				bash ./Anaconda3-5.2.0-Linux-x86_64.sh
				source ~/.bashrc
			fi
		fi
		if ! [ -f ~/.bashrc ]; then
			touch ~/.bashrc
			if ! [ -x "$(command -v wget)" ]; then
				echo 'Kindly install "Wget" and run this script again, thanks!!' >&2
				exit 1
			fi
			if [ -x "$(command -v wget)" ]; then
				wget https://repo.anaconda.com/archive/Anaconda3-5.2.0-Linux-x86_64.sh
				bash ./Anaconda3-5.2.0-Linux-x86_64.sh
				source ~/.bashrc
			fi
		fi
	fi
fi

if [ "$(uname)" = "Linux" ]; then
	if [ "$(uname -m)" != "x86_64" ]; then
		if [ -f ~/.bashrc ]; then
			if ! [ -x "$(command -v wget)" ]; then
				echo 'Kindly install "Wget" and run this script again, thanks!!' >&2
				exit 1
			fi
			if [ -x "$(command -v wget)" ]; then
				wget https://repo.anaconda.com/archive/Anaconda3-5.2.0-Linux-x86.sh
				bash ./Anaconda3-5.2.0-Linux-x86.sh
				source ~/.bashrc
			fi
		fi
		if ! [ -f ~/.bashrc ]; then
			touch ~/.bashrc
			if ! [ -x "$(command -v wget)" ]; then
				echo 'Kindly install "Wget" and run this script again, thanks!!' >&2
				exit 1
			fi
			if [ -x "$(command -v wget)" ]; then
				wget https://repo.anaconda.com/archive/Anaconda3-5.2.0-Linux-x86.sh
				bash ./Anaconda3-5.2.0-Linux-x86.sh
				source ~/.bashrc
			fi
		fi
	fi
fi

if [ -x "$(command -v conda)" ]; then
	tset
	reset
	echo 'Conda Installed Successfully!!' >&2
fi

if [ -x "$(command -v conda)" ]; then
	conda install -c biocore hmmer
	conda install -c conda-forge/label/gcc7 biopython
	conda install -c etetoolkit ete3 ete_toolchain
	ete3 build check
fi

if [ -x "$(command -v conda)" ]; then
	if [ -x "$(command -v jackhmmer)" ]; then
		if [ -x "$(command -v python -c "import Bio")" ]; then
			if [ -x "$(command -v ete3)" ]; then
				echo 'FlaGs Environment Installed Successfully!!' >&2
				exit 1
			fi
			echo 'FlaGs Environment Not Installed Properly!!' >&2
			exit 1
		fi
		echo 'FlaGs Environment Not Installed Properly!!' >&2
		exit 1
	fi
	echo 'FlaGs Environment Not Installed Properly!!' >&2
	exit 1
fi
echo 'FlaGs Environment Not Installed Properly!!' >&2
exit 1
