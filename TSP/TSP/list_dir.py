import sys
import os
import re

def dir_exist(dir):
    if(os.path.isdir(dir)):
        return dir
    else:
        print("You must insert a dir as argument or in the script input")
        sys.exit()


if len(sys.argv) <= 3:
    if len(sys.argv) >=2:
        dir = dir_exist(sys.argv[1])

        if len(sys.argv) ==3:
            max_n = int(sys.argv[2])

            if(max_n==-1):
                max_n=90000

    else:
        dir = input()
        dir = dir_exist(dir)
        max_n = int(input())

        if(max_n==-1):
            max_n=90000

else:
    print("Too many arguments")
    sys.exit()

print(dir)
# print([x for x in os.listdir(sys.argv[1]) if x.endswith('.tsp')])

file = open("instances.txt", "w")
[file.write(dir+"/"+x+"\n")  for x in os.listdir(dir) if (x.endswith('.tsp') and int(re.findall(r'\d+',x)[0])<=max_n)]
#[file.write(x+"\n") for x in tsp_files]
file.close()
