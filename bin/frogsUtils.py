#
# Copyright (C) 2014 INRA
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__author__ = 'Frederic Escudie - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '0.2.0'
__email__ = 'frogs@toulouse.inra.fr'
__status__ = 'prod'

import os
import sys
import time
import subprocess
from subprocess import Popen, PIPE


def which(exec_name):
    """
    @summary: Returns the software absolute path.
    @param exec_name: [str] The software file (example : blastn or classifier.jar).
    @return: [str] The software absolute path.
    """
    path_locations = os.getenv("PATH").split(os.pathsep) + sys.path
    exec_path = None
    for current_location in path_locations:
        if exec_path is None and os.path.isfile(os.path.join(current_location, exec_name)):
            exec_path = os.path.abspath( os.path.join(current_location, exec_name) )
    if exec_path is None:
        raise Exception( "The software '" + exec_name + "' cannot be retrieved in path." )
    return exec_path


def prevent_shell_injections(argparse_namespace, excluded_args=None):
    """
    @summary: Raises an exception if one parameter contains a backquote or a semi-colon.
    @param argparse_namespace: [Namespase] The result of parser.parse_args().
    @param excluded_args: [list] List of unchecked parameters.
    """
    exceptions = list() if excluded_args is None else excluded_args
    for param_name in argparse_namespace.__dict__.keys():
        if not param_name in exceptions:
            param_val = getattr(argparse_namespace, param_name)
            if issubclass(param_val.__class__, list):
                new_param_val = list()
                for val in param_val:
                    if ';' in val.encode('utf8') or '`' in val.encode('utf8') or '|' in val.encode('utf8'):
                        raise Exception( "';' and '`' are unauthorized characters." ) 
            elif param_val is not None and issubclass(param_val.__class__, str):
                if ';' in param_val.encode('utf8') or '`' in param_val.encode('utf8') or '|' in param_val.encode('utf8'):
                    raise Exception( "';' and '`' are unauthorized characters." )


class Cmd:
    """
    @summary : Command wrapper.
    """
    def __init__(self, program, description, exec_parameters, version_parameters=None):
        """
        @param exec_parameters: [str] The parameters to execute the program. Two possibles syntaxes.
                                If the parameter contains the string '##PROGRAM##', this tag will be replaced by the program parameter before submit.
                                Otherwise the parameters will be added after the program in command line.
        @param version_parameters: [str] The parameters to get the program version. Two possibles syntaxes.
                                   If the parameter contains the string '##PROGRAM##', this tag will be replaced by the program parameter before submit.
                                   Otherwise the parameters will be added after the program in command line.
        """
        self.program = program
        self.description = description
        self.exec_parameters = exec_parameters
        self.version_parameters = version_parameters

    def get_cmd(self):
        """
        @summary : Returns the command line.
        @return : [str] The command line.
        """
        cmd = None
        if '##PROGRAM##' in self.exec_parameters:
            cmd = self.exec_parameters.replace('##PROGRAM##', self.program)
        else:
            cmd = self.program + ' ' + self.exec_parameters
        return cmd

    def get_version(self, location='stderr'):
        """
        @summary : Returns the program version number.
        @param location : [str] If the version command returns the version number on 'stdout' or on 'stderr'.
        @return : [str] version number if this is possible, otherwise this method return 'unknown'.
        """
        if self.version_parameters is None:
            return "unknown"
        else:
            try:
                cmd = self.program + ' ' + self.version_parameters
                if '##PROGRAM##' in self.exec_parameters:
                    cmd = self.version_parameters.replace('##PROGRAM##', self.program)
                p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
                stdout, stderr = p.communicate()
                if location == 'stderr':
                    return stderr.strip()
                else:
                    return stdout.strip()
            except:
                raise Exception( "Version cannot be retrieve for the software '" + self.program + "'." )

    def parser(self, log_file):
        """
        @summary : Parse the command results to add information in log_file.
        @log_file : [str] Path to the sample process log file.
        """
        pass

    def submit(self, log_file=None):
        """
        @summary : Launch command, trace this action in log and parse results.
        @log_file : [str] Path to the sample process log file.
        """
        # Log
        if log_file is not None:
            FH_log = Logger( log_file )
            FH_log.write( '# ' + self.description + ' (' + os.path.basename(self.program) + ' version : ' + self.get_version() + ')\n' )
            FH_log.write( 'Command:\n\t' + self.get_cmd() + '\n' )
            FH_log.write( 'Execution:\n\tstart: ' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + '\n' )
            FH_log.close()
        # Process
        subprocess.check_output( self.get_cmd(), shell=True )
        # Log
        if log_file is not None:
            FH_log = Logger( log_file )
            FH_log.write( '\tend:   ' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + '\n' )
            FH_log.close()
            # Post-process results
            self.parser(log_file)


class Logger:
    """
    @summary: Log file handler.
    """

    def __init__(self, filepath=None):
        """
        @param filepath: [str] The log filepath. [default : STDOUT]
        """
        self.filepath = filepath
        if self.filepath is not None and self.filepath is not sys.stdout:
            self.file_handle = open( self.filepath, "a" )
        else:
            self.file_handle = sys.stdout

    def __del__(self):
        """
        @summary: Closed file handler when the logger is detroyed.
        """
        self.close()

    def close(self):
        """
        @summary: Closed file handler.
        """
        if self.filepath is not None and self.filepath is not sys.stdout:
            if self.file_handle is not None:
                self.file_handle.close()
                self.file_handle = None

    def write(self, msg):
        """
        @summary: Writes msg on file.
        @param msg: [str] The message to write.
        """
        self.file_handle.write( msg )

    @staticmethod
    def static_write(filepath, msg):
        """
        @summary: Writes msg on file.
        @param filepath: [str] The log filepath. [default : STDOUT]
        @param msg: [str] The message to write.
        """
        if filepath is not None and filepath is not sys.stdout:
            FH_log = open( filepath, "a" )
            FH_log.write( msg )
            FH_log.close()
        else:
            sys.stdout.write( msg )


class TmpFiles:
    """
    @summary: Manager for temporary files.
    @note:
        tmpFiles = TmpFiles(out_dir)
        try:
            ...
            tmp_seq = tmpFiles.add( "toto.fasta" )
            ...
            tmp_log = tmpFiles.add( "log.txt" )
            ...
        finaly:
            tmpFiles.deleteAll()
    """
    def __init__(self, tmp_dir, prefix=None):
        """
        @param tmp_dir: [str] The temporary directory path.
        @param prefix: [str] The prefix added to each temporary file [default: <TIMESTAMP>_<PID>].
        """
        if prefix is None:
            prefix = str(time.time()) + "_" + str(os.getpid())
        self.files = list()
        self.dirs = list()
        self.tmp_dir = tmp_dir
        self.prefix = prefix

    def add(self, filename, prefix=None, dir=None):
        """
        @summary: Add a temporary file.
        @param filename: The filename without prefix.
        @param prefix: The prefix added [default: TmpFiles.prefix].
        @param dir: The directory path [default: TmpFiles.tmp_dir].
        @return: [str] The filepath.
        """
        # Default
        if prefix is None:
            prefix = self.prefix
        if dir is None:
            dir = self.tmp_dir
        # Process
        filepath = os.path.join(dir, prefix + "_" + filename)
        self.files.append(filepath)
        return filepath

    def add_dir(self, dirname, prefix=None, dir=None):
        """
        @summary: Add a temporary dir.
        @param filename: The dirname without prefix.
        @param prefix: The prefix added [default: TmpFiles.prefix].
        @param dir: The directory path [default: TmpFiles.tmp_dir].
        @return: [str] The filepath.
        """
        # Default
        if prefix is None:
            prefix = self.prefix
        if dir is None:
            dir = self.tmp_dir
        # Process
        dirpath = os.path.join(dir, prefix + "_" + dirname)
        self.dirs.append(dirpath)
        return dirpath

    def delete(self, filepath):
        """
        @summary: Deletes the specified temporary file.
        @param filepath: [str] The file path to delete.
        """
        self.files.remove(filepath)
        if os.path.exists(filepath): os.remove(filepath)

    def delete_dir(self, dirpath):
        """
        @summary: Deletes the specified temporary dir.
        @param filepath: [str] The file path to delete.
        """
        if dirpath in self.dirs: self.dirs.remove(dirpath)
            
        if os.path.exists(dirpath):
            for root, dirnames,filenames in os.walk(dirpath):
                for f in filenames:
                    if f in self.files: self.files.remove(os.path.join(dirpath,f))
                    if os.path.exists(os.path.join(dirpath,f)): os.remove(os.path.join(dirpath,f))
                for d in dirnames:
                    if d in self.dirs: self.dirs.remove(os.path.join(dirpath,d))
                    if os.path.exists(os.path.join(dirpath,d)): self.delete_dir(os.path.join(dirpath,d))
            os.rmdir(dirpath)

    def deleteAll(self):
        """
        @summary: Deletes all temporary files.
        """
        all_tmp_files = [tmp_file for tmp_file in self.files]
        for tmp_file in all_tmp_files:
            self.delete(tmp_file)
        
        all_tmp_dirs=[tmp_dir for tmp_dir in self.dirs]
        for tmp_dir in all_tmp_dirs:
            self.delete_dir(tmp_dir)