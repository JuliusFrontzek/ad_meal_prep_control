from setuptools import setup
import os


setup(
    name="ad_meal_prep_control",
    version="0.1.0",
    description="Multi-stage Model Predictive Control for biogas plant \
        (anaerobic digestion a.k.a. AD) to optimally control substrate feed (meal prep).",
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
)

os.mkdir("./results/")
