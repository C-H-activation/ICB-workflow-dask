# dags/project/package/bash.py

try:
    import os
except ImportError:
    raise ImportError(
        "Module os is required. Please install module os."
    )
try:
    import shutil
except ImportError:
    raise ImportError(
        "Module shutil is required. Please install module shutil."
    )

from project.package.utility import existDirectory

def mkdir_p(path: str):
    """Simulate bash command 'mkdir -p'."""
    try:
        os.makedirs(path)
    except OSError as error:
        if os.path.isdir(path):
            pass
        else:
            raise error


def createDirectory(name: str, path=os.getcwd()) -> str:
    """Create directory for mentioned name."""

    path = os.path.join(path, name)

    # create temporary directory
    if not os.path.isdir(path):
        mkdir_p(path)

    return path


def copyFile(source: str, destination: str):
    """Copy file from source (path) to destination (path)."""
    shutil.copy2(source, destination)


def changeDirectory(destination: str):
    """Change directory to destination."""
    if existDirectory(destination):
        os.chdir(destination)
