import sys
import os

def dir_exist(dir):
    if(os.path.isdir(dir)):
        return dir
    else:
        print("You must insert a dir as argument or in the script input")
        sys.exit()


if len(sys.argv) == 2:
    dir = dir_exist(sys.argv[1])
elif len(sys.argv) > 2:
    print("Too many arguments")
    sys.exit()
else:
    dir = input();
    dir = dir_exist(dir)

print(dir)
# print([x for x in os.listdir(sys.argv[1]) if x.endswith('.tsp')])

file = open("instances.txt", "w")
[file.write(dir+"/"+x+"\n")  for x in os.listdir(dir) if x.endswith('.tsp')]
#[file.write(x+"\n") for x in tsp_files]
file.close()
