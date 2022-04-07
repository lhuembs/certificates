In order to run the program you need a license of GUROBI and the miplibing library from
https://github.com/thserra/MIPLIBing.
In a second step clone the conda environment in environment.yml and activate it.
The documentation for conda can be found [here](https://docs.conda.io/en/latest/)
and you can find the commands for creating an environment from the the .yml file [here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#create-env-from-file)

Create a directory "Certificates". Here the certificates will be saved.

Then the program can be used with the command
`python compute_certificate.py flugpl`
in order to compute a certificate to the MIPLIB problem 'flugpl'.
When computing other problems from the MIPLIB just replace the name.
If you want to compute a certificate for a problem where you have the .lp or .mps file locally use the the options `--local` or `--mps` respectively.
The program also provides a `--help` option for a more detailed description.

You can find the corresponding preprint to the program [here](https://opus4.kobv.de/opus4-trr154/frontdoor/index/index/docId/476/).
