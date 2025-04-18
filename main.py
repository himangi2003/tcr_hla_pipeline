# main.py
from config import *
from train import run_training
from test import run_testing
from filtering_data import run_filtering

if __name__ == "__main__":
    edit_type = "edit1"
    min_num_ptids = 4
    only_novel = True
    os.makedirs("outputs", exist_ok=True)

    print("Filtering and Preprocessing data....")
    run_filtering(min_num_ptids, only_novel)

    print("Running training...")
    run_training(edit_type, min_num_ptids, only_novel)

    print("Running testing...")
    run_testing(edit_type)

