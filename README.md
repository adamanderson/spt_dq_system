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
  openssl rsa -passin pass:x -in server.pass.key -out server.key
  rm server.pass.key
  openssl req -new -key server.key -out server.csr
  openssl x509 -req -sha256 -days 365 -in server.csr -signkey server.key -out server.crt
  ```

* bcrypt hash
  * A file containing a bcrypt hash of the website password is needed. Generate the hash by placing the following code in db_server.js
  ```javascript
  var salt = bcrypt.genSaltSync(10);
  var hash = bcrypt.hashSync("B4c0/\/", salt);
  console.log(hash);
  exit();
  ```
  * Then copy the hash from the terminal into a file called "hash" and place in the server's directory.

# Running

To launch the server on Scott or Amundsen:
```bash
node db_server.js
```
This will start a server at 127.0.0.1 on port 3000. Then point your browser at http://127.0.0.1:3000/index.html.

# More Information
For a more thorough description of how to use and extend the website, please visit the wiki page at https://pole.uchicago.edu/spt3g/index.php/SPT_Data_Quality
