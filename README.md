# nanoPolish
Polish a nanopore assembly using Racon and Medaka



The hd5 files that are used as models are stored in gits large file storage `lfs`, if you do a git clone those files will just be text pointers unless you have git lfs installed.  If you run into a `Unable to open file (file signature not found)` error it is because the hd5 file is not an actual model file but a pointer to git's lfs storage.

```
git clone git@github.com:nanoporetech/medaka.git
git lfs install
git lfs pull
```
