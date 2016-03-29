#!/usr/bin/env python
"""
Usage:
    WP-BACKUP.PY WPFOLDER

    WPFOLDER is the root folder of your Wordpress installation, that is, it
    contains the file wp-config.php.

The script archives all important wordpress data and configurations in the
following archives where YYYY-MM-DD represents the current date. Archive files
are created in the current working directory.

YYYY-MM-DD-config.zip
    Stores the wp-config.php file which holds things like database connection
    parameters and secure hashes.
YYYY-MM-DD-content.zip
    Stores the wp-content folder which holds translations, themes, plugins, and
    that kind of stuff. The uploads folder is excluded.
YYYY-MM-DD-uploads.zip
    Stores the wp-content/uploads folder separately from the other wp-content.
YYYY-MM-DD-sql.zip
    Stores a dump of the database. This contains things like all your posts,
    the comments, Wordpress and plugin configurations, user accounts, and so on.


Note:
    Unfortunately the implementation of ZipFile in Python does not have any
    means of password protecting the files. Given the sensitive nature of the
    information stored, it is advised to run this script in a folder where only
    the current user has reading permissions.

Update (27-08-2013):
	Archives now have only read permission to the user.
Update (20-02-2014):
	Should now be compatible with Python 3.
"""
# Author:   Paul Koppen
# Website:  http://paulkoppen.com/
# Email:    wp-backup () paulkoppen.com
# Date:     15 November 2011

#from __future__ import with_statement
from contextlib import closing
from datetime   import date
from pexpect    import spawn, EOF, TIMEOUT
from stat       import S_IRUSR
from tempfile   import mkstemp
from zipfile    import ZipFile, ZIP_DEFLATED
import re
import os



def backupwordpress(wpfolder):
	''' backup wordpress files
	the following files and folders are stored:
		- wp-content folder (separate file for uploads subfolder)
		- wp-config.php with connection parameters and hashes
		- the complete database
	in most circumstances all other files are part of the downloaded wordpress
	installation and therefore require no backup. if they do you're probably
	better off simply zipping the entire wordpress folder.
	'''
	# --- wp-config.php
	configfile  = os.path.join(wpfolder, 'wp-config.php')
	archivename = zipname('config')
	with closing(ZipFile(archivename, 'w', ZIP_DEFLATED)) as archive:
		addtozip(configfile, archive)
	os.chmod(archivename, S_IRUSR)
	
	# --- wp-content/		(excluding uploads)
	folder      = os.path.join(wpfolder, 'wp-content')
	contents    = (os.path.join(folder,fname)
				for fname in os.listdir(folder) if fname != 'uploads')
	archivename = zipname('content')
	with closing(ZipFile(archivename, 'w', ZIP_DEFLATED)) as archive:
		for path in contents:
			addtozip(path, archive, 'wp-content')
	os.chmod(archivename, S_IRUSR)
	
	# --- wp-content/uploads/
	folder      = os.path.join(wpfolder, 'wp-content', 'uploads')
	archivename = zipname('uploads')
	with closing(ZipFile(archivename, 'w', ZIP_DEFLATED)) as archive:
		addtozip(folder, archive, 'wp-content')
	os.chmod(archivename, S_IRUSR)
	
	# --- SQL
	archivename = zipname('sql')
	sqlfname    = archivename[:-4]
	tmpfname    = dumpsql(*readwpconfig(wpfolder))
	os.rename(tmpfname, sqlfname)
	with closing(ZipFile(archivename, 'w', ZIP_DEFLATED)) as archive:
		addtozip(sqlfname, archive)
	os.unlink(sqlfname)
	os.chmod(archivename, S_IRUSR)



def zipname(basename):
	''' create a suitable name for a zipfile from a basename, e.g. foldername
	'''
	filename = os.path.basename(basename) or \
				os.path.basename(os.path.dirname(basename))
	filename = filename[3:] if filename.startswith('wp-') else filename
	today    = date.today()
	return '%04d-%02d-%02d.%s.zip' % (
				today.year, today.month, today.day, filename)

def addtozip(path, archive, relfolder=''):
	''' add a file or folder to a ZipFile archive
	optionally prepend a relative folder, so the contents are "put in that folder"
	@note ignores empty directories
	'''
	abspath      = os.path.abspath(path)
	parent,fname = os.path.split(abspath)
	if os.path.isdir(abspath):
		offset = len(parent + os.path.sep)
		for root, dirs, files in os.walk(abspath):
			for fname in files:
				absname = os.path.join(root, fname)
				relname = os.path.join(relfolder, absname[offset:])
				archive.write(absname, relname)
	else:
		archive.write(abspath, os.path.join(relfolder, fname))

def dumpsql(database, username, password, hostname):
	''' dump database to temporary file in cwd
	return filename
	'''
	# note: do not send password as command-line argument because it will be
	#       readable to anyone on the machine using `ps aux | grep mysqldump`.
	command          = 'mysqldump'
	arguments        = ['-h',hostname, '-u',username, '-p', database]
	(tmpfd,tmpfname) = mkstemp(dir=os.path.abspath('./'), text=True)
	tmpfile          = os.fdopen(tmpfd, 'w')
	childproc        = spawn(command, arguments)
	childproc.expect('Enter password:')
	childproc.sendline(password)
	childproc.logfile = tmpfile
	while childproc.expect([EOF, TIMEOUT]):
		pass
	childproc.close()
	tmpfile.close()
	return tmpfname

def readwpconfig(wpfolder):
	''' read defined database connection parameters wp-config.php
	return tuple of values in the order expected by dumpsql
	'''
	configfile = os.path.join(wpfolder, 'wp-config.php')
	command    = 'php'
	arguments  = ['-r',('include(\'%s\');var_dump(DB_NAME,DB_USER,DB_PAS' +\
				'SWORD,DB_HOST);')%configfile]
	# e.g. 'string(10) "abcdefghij"\r\n'
	re_var     = re.compile(r'string\(\d+\) "(.*)"')
	childproc  = spawn(command, arguments)
	return tuple(re.match(re_var,line).group(1) for line in childproc)



if __name__ == '__main__':
	import sys
	if len(sys.argv) < 2 or sys.argv[1] in ('/?','/h','-h','--help'):
		print(__doc__)
		sys.exit(1)
	basedir = sys.argv[1]
	if basedir.endswith(os.path.sep):
		basedir = basedir[:-len(os.path.sep)]
	backupwordpress(basedir)


