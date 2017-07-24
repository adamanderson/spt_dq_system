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

* To install node.js module dependencies, run on the server:
  ```bash
  cd spt_dq_system
  npm install
  ```
* openssl
  * Create an openssl security key in order to use https functionality.
  ```bash
  openssl genrsa -des3 -passout pass:x -out server.pass.key 2048
  openssl rsa -passin pass:x -in server.pass.key -out key.pem
  rm server.pass.key
  openssl req -new -key key.pem -out server.csr
  openssl x509 -req -sha256 -days 365 -in server.csr -signkey key.pem -out cert.pem
  rm server.csr
  ```

* bcrypt hash
  * A file containing a bcrypt hash of the website password is needed. Generate the hash by placing the following code in db_server.js
  ```javascript
  var salt = bcrypt.genSaltSync(10);
  var hash = bcrypt.hashSync("your_password", salt);
  console.log(hash);
  exit();
  ```
  * Then copy the hash from the terminal into a file called "hash" and place in the server's directory.

# Running

To launch the server on Scott or Amundsen:
```bash
node db_server.js
```
This will start a server at https://127.0.0.1:3000 if on Scott or https://scott.grid.uchicago.edu:3000 if on a remote machine.

# More Information
For a more thorough description of how to use and extend the website, please visit the wiki page at https://pole.uchicago.edu/spt3g/index.php/SPT_Data_Quality
