# main.py
from config import *
from train import run_training
from test import run_testing

if __name__ == "__main__":
    edit_type = "edit1"
    min_num_ptids = 2
    only_novel = True

    print("Running training...")
    run_training(edit_type, min_num_ptids, only_novel)

    print("Running testing...")
    run_testing(edit_type)

