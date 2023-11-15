import argparse
import os
import logging

def main():
    description = "Script to download the ProMGE database"
    softwareName = "promgedl"
    parser = argparse.ArgumentParser(
        prog=softwareName,
        description=description,
        usage=f"{softwareName} [option] ...",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

if __name__ == "__main__":
    main()