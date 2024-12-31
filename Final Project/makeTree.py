import os

def write_tree_structure(start_path, output_file, indent=""):
    """Recursively writes the folder tree structure starting from start_path to a file."""
    try:
        # Get a list of items in the directory
        items = os.listdir(start_path)
        for i, item in enumerate(sorted(items)):
            item_path = os.path.join(start_path, item)
            # Check if the item is the last one
            is_last = i == len(items) - 1
            # Write with appropriate tree characters
            output_file.write(f"{indent}{'└── ' if is_last else '├── '}{item}\n")
            # Recurse into directories
            if os.path.isdir(item_path):
                write_tree_structure(item_path, output_file, indent + ("    " if is_last else "│   "))
    except PermissionError:
        output_file.write(f"{indent}[Permission Denied]\n")

if __name__ == "__main__":
    # Automatically set the root directory to the current working directory
    root_directory = os.getcwd()
    output_file_path = "tree.txt"

    with open(output_file_path, "w") as output_file:
        output_file.write(f"Root Directory: {root_directory}\n")
        write_tree_structure(root_directory, output_file)

    print(f"Tree structure has been written to {output_file_path}")
