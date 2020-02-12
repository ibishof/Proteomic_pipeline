import sys
import suds

_script = 'GO_Elite.py'
_appName = "GO_Elite"
_appVersion = '1.2.4'
_appDescription = "GO-Elite (http://genmapp.org/go_elite) is a software tool designed to identify a minimal non-redundant set "
_appDescription +="of Gene Ontology (GO) biological terms or pathways to describe a particular set of genes. In addition, "
_appDescription +="alternate ontologies (e.g., Disease Ontology), gene sets and metabolomics data can also be used as input."
_authorName = 'Nathan Salomonis'
_authorEmail = 'nsalomonis@gmail.com'
_authorURL = 'http://genmapp.org/go_elite'
_appIcon = "goelite.ico"

excludes = ["matplotlib", "wxPython"] #"numpy","scipy",
includes = ["suds"]
""" By default, suds will be installed in site-packages as a .egg file (zip compressed). Make a duplicate, change to .zip and extract
here to allow it to be recognized by py2exe (must be a directory) """

#data_files=matplotlib.get_py2exe_datafiles()

matplot_exclude = ['MSVCP90.dll']
scipy_exclude = ['libiomp5md.dll','libifcoremd.dll','libmmd.dll']

""" xml.sax.drivers2.drv_pyexpat is an XML parser needed by suds that py2app fails to include. Identified by looking at the line: parser_list+self.parsers in
/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/PyXML-0.8.4-py2.7-macosx-10.6-intel.egg/_xmlplus/sax/saxexts.py
check the py2app print out to see where this file is in the future """

if sys.platform.startswith("darwin"):
        ### example command: python setup.py py2app
        includes+= ["xml.sax.drivers2.drv_pyexpat"]
        from distutils.core import setup
        import py2app
        options = {"py2app":
                    {"excludes": excludes,
                     "includes": includes,
                     #argv_emulation = True,
                    "iconfile": "goelite.icns"}
        }
        setup(name=_appName,
                        app=[_script],
                        version=_appVersion,
                        description=_appDescription,
                        author=_authorName,
                        author_email=_authorEmail,
                        url=_authorURL,
                        options=options,
                        #data_files=data_files,
                        setup_requires=["py2app"]
        )

if sys.platform.startswith("win"):
        ### example command: python setup.py py2exe
        from distutils.core import setup
        import py2exe
        import suds
        windows=[{"script":_script,"icon_resources":[(1,_appIcon)]}]
        options={'py2exe':
                        {
                        "includes": 'suds',
                        #"includes": 'matplotlib',
                        #"includes": 'mpl_toolkits',
                        "dll_excludes": matplot_exclude+scipy_exclude,
                        }}
        setup(
                        windows = windows,
                        options = options,
                        version=_appVersion,
                        description=_appDescription,
                        author=_authorName,
                        author_email=_authorEmail,
                        url=_authorURL,
                        #data_files=data_files,
        )

if sys.platform.startswith("2linux"):
	# bb_setup.py
	from bbfreeze import Freezer
	 
	f = Freezer(distdir="bb-binary")
	f.addScript("GO_Elite.py")
	f()

if sys.platform.startswith("linux"):
        ### example command: python setup.py build
        from cx_Freeze import setup, Executable
        ### use to get rid of library.zip and move into the executable, along with appendScriptToLibrary and appendScriptToExe
        #buildOptions = dict(create_shared_zip = False) 

	setup(
		name = _appName,
		version=_appVersion,
		description=_appDescription,
		author=_authorName,
		author_email=_authorEmail,
		url=_authorURL,
		#options = dict(build_exe = buildOptions),
		executables = [Executable(_script,
				#appendScriptToExe=True,
				#appendScriptToLibrary=False,
				#icon='goelite.ico',
				compress=True)],
	)

