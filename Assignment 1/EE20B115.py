import sys                           # importing library

n = len(sys.argv)                    # counting the number of arguments provided

try:                                 # start of the try block

 if (n == 1):
    print("\nThe name of the netlist file is not provided.\n")            # if only name of python file is provided, then this message
 elif(n == 2):
    print("\nThe name of the netist file provided "+ sys.argv[1]+".\n")   # using 2nd argument as name of netlist file

 filename = sys.argv[1]          

 f = open(filename)                        # opening file in python with f as the file object

 lines = f.readlines()                     # using readlines method of f (file object) to return each line of file as an entry in a list named lines            
      
 i = 0
 while (lines[i]!=".circuit\n"):           # finding the position of .circuit in list
  i = i+1

 j = 0
 while ((lines[j]!=".end") and (lines[j]!=".end\n")):    # finding the position of .end in list
  j = j+1

 print("The netlist data in reverse order is:\n")

 for k in range(j-1,i,-1):                        # looping in reverse order from higher to lower value
  m = lines[k].split('#')[0].split()              # Kth element of lines list is a string which is one of the lines of netlist,  
                                                  # first splitting the string by '#' and taking the part before 1st occurence of'#' using [0]
                                                  # then splitting that further by whitepsaces 

  m.reverse()                                            # reversing the elements of the list
  n = " ".join(m)                                        #  converting the list into a string with elements seperated by whitespace
  print(n)                                       

 f.close()                                               # closing the file after printing the reversed order

except IOError:                                          # for Input/Output operation error, given message is displayed
 print("Incorrect file name or file path.")
 sys.exit()

except IndexError:                                       # for Indexing error, given message is displayed
  print("Invalid netlist format.") 
  sys.exit()
 
