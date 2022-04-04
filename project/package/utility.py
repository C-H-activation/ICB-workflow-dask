# dags/project/package/utility.py

try:
    import os
except ImportError:
    raise ImportError(
        "Module os is required. Please install module os."		
    )


def existFile(path: str) -> bool:
    """Check if file with 'path' exists."""
    return os.path.isfile(path)


def existDirectory(path: str) -> bool:
    """Check if directory with 'path' exists."""
    return os.path.isdir(path)

def is_exe(path):
	"""Check if program is availanle and executable."""
	return os.path.isfile(path) and os.access(path, os.X_OK)

def program(name: str) -> bool:
    """Check if a program is in path and executable."""
    
    fpath, fname = os.path.split(name)
    
    if fpath:
        if is_exe(name):
            return True
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            executable = os.path.join(path, name)
            if is_exe(executable):
                return True
    return False
