Data quality monitoring system for SPT3G.

# Setup
The following software is required on the server side:
* node.js
  * Get it from the node.js website: https://nodejs.org/en/download/
  * Or better yet, from your Linux package manager:
  ```bash
  apt-get install nodejs
  ```
  * ... or homebrew if you traffick in Mac OS:
  ```bash
  brew install node
  ```

To install node.js module dependencies, run on the server:
```bash
cd spt_dq_system
npm install
```

# Running
First build the database on scott.uchicago.edu:
```bash
python build_database.py
```

This will create a database file `obsfiles.db`. Copy this to wherever you want to run the webserver (currently this must be your local machine). To launch the server on your local machine:
```bash
node db_server.js
```
This will start a server at 127.0.0.1 on port 3000. Then point your browser at http://127.0.0.1:3000/index.html.

# Adding More Plots
To add new plot types to the webserver, simply create a new python script in the plot directory subject to the following constraints:
* The name of the plot type is the title of the file (minus the '.py')
* The script contains a function with the same name as the file (minus the '.py')
* That function returns a matplotlib figure
* That function takes one argument that is a dictionary containing the source, observation id, path, and plot type

Look in the plot directory for examples of the format.
