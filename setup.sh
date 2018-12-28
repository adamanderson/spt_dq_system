# check arguments
if [ "$1" != "pole" ] && [ "$2" != "north" ];
then
	echo "Did not specify site argument of 'pole' or 'north'. "\
	"Refusing to run."
	exit 1
fi

# get node.js dependencies (server-side)
npm install

# get local copy of js dependencies (client-side)
cd js
wget https://ajax.googleapis.com/ajax/libs/jqueryui/1.12.1/themes/smoothness/jquery-ui.css
mv jquery-ui.css ../css
wget https://code.jquery.com/jquery-1.12.4.js
wget https://code.jquery.com/ui/1.12.1/jquery-ui.js
wget https://momentjs.com/downloads/moment.min.js
wget https://cdnjs.cloudflare.com/ajax/libs/tabulator/3.1.5/css/tabulator.min.css
mv tabulator.min.css ../css
wget https://cdnjs.cloudflare.com/ajax/libs/tabulator/3.1.5/js/tabulator.min.js
wget https://cdnjs.cloudflare.com/ajax/libs/gsap/1.18.2/TweenMax.min.js
cd ..

# clone spt3g_software
if [ "$1" == "pole" ];
then
	git clone sptviewer@anal:sptRepos/spt3g_software
elif [ "$1" == "north" ];
then
	git clone https://github.com/SouthPoleTelescope/spt3g_software.git
else
	echo "Did not specify site argument of 'pole' or 'north'. "\
	"Cannot clone spt3g_software."
	exit 1
fi

# build spt3g_software
mkdir spt3g_software/build
cd spt3g_software/build
cmake ..
make -j4
cd ..
