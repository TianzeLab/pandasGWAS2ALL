import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pandasgwas",
    version="2.0.0",
    author="Cao Tianze",
    author_email="hnrcao@qq.com",
    description="A Python package for easy retrieval of GWAS Catalog data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/TianzeLab/pandasGWAS2",
    project_urls={
        "Bug Tracker": "https://github.com/TianzeLab/pandasGWAS2/issues",
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.12",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    # package_dir={"": "pandasgwas"},
    packages=['pandasgwas'],
    python_requires=">=3.12",
    install_requires=['pandas>=3.0.3', 'requests>=2.34.2', 'progressbar2>=4.5.0', 'cachetools>=7.1.4'],
    license="MIT",
    keywords=['gwas', 'genomics', 'snp', 'bioinformatics','pandas'],
    package_data={
        "": ["*.csv","*.txt"]
    }
)
