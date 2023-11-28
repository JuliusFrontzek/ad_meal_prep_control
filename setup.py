from setuptools import setup

setup(
    name="ad_meal_prep_control",
    version="0.1.0",
    description="Multi-stage Model Predictive Control algorithm and simulation \
        for a biogas plant (anaerobic digestion) with an optional external gas storage. \
        Considers uncertainties in macronutrients of the fed substrates.",
    long_description=open("README.md").read(),
    author="Julius Frontzek",
    author_email="ju-frontzek@gmx.de",
    packages=["ad_meal_prep_control"],
    install_requires=[
        "numpy",
        "matplotlib",
        "do-mpc",
        "uncertainties",
        "pygame",
    ],
    data_files=[("./results", [])],
)
