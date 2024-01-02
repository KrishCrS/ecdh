# Elliptic Curve Diffie-Hellman
Just run the program. All things are set. Also if you wish you can change the secret values  `a`  of Alice's or  `b`  of Bob's. We just check whether a secret private key is generated or not. If   `a*(b*P) == b*(a*P)`  then we are succeed and secure private key  `abP` is generated and they can use it as for further coomunications.  

## Run Makefile :
1. To make : `make`
2. To make and run : `make run`
3. To clean : `make clean`

### OR

## Run the following commands :
1.  `gcc ecdh.c functions.c -o output`
2.  `./output`
