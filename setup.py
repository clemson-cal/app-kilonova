import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="knc_tools",
    version="0.1.0",
    author="J. Zrake",
    author_email="jzrake@clemson.edu",
    description="Python module to work with data output from the kilonova code",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/clemson-cal/app-kilonova",
    packages=setuptools.find_packages(),
    entry_points={
        'console_scripts': [
            'knc-plot = knc_tools.plot:main',
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=['numpy', 'matplotlib'],
    python_requires='>=3.6',
)
