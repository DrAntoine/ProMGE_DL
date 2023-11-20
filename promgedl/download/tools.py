import os



def check_directory_structure(structure):
    folderStatus = []
    for s in structure:
        folderStatus.append(os.path.isdir(s))
    if len(set(folderStatus))==2:
        return "Partial"
    else:
        if list(set(folderStatus))[0]==True:
            return "Complete"
        else:
            return "Absent"