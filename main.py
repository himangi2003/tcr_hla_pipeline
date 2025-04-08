# main.py
from config import *
from train import run_training
from test import run_testing

if __name__ == "__main__":
    print("Running training...")
    run_training()

    print("Running testing...")
    run_testing()

