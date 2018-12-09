from setuptools import setup
from setuptools import find_packages

def setup_package():
  install_requires = ['pandas', 'scipy', 'numpy']
  # console_scripts = [ 'biofam=biofam.build_model.init_model:entry_point'],
  metadata = dict(
      name = 'svca',
      version = '0.0.1',
      author = 'Damien Arnol',
      packages = find_packages(),
      install_requires = install_requires,
    )

  setup(**metadata)

if __name__ == '__main__':
  setup_package()
