hi there!!!
and hello!

---

einsteinium$ git clone https://github.com/Macaulay2/Workshop-2015-Boise.git
Cloning into 'Workshop-2015-Boise'...
remote: Counting objects: 18, done.        
remote: Compressing objects:  50% (1/2)           
remote: Compressing objects: 100% (2/2)           
remote: Compressing objects: 100% (2/2), done.        
Unpacking objects:   5% (1/18)   
Unpacking objects:  11% (2/18)   
Unpacking objects:  16% (3/18)   
Unpacking objects:  22% (4/18)   
remote: Total 18 (delta 0), reused 0 (delta 0), pack-reused 14        
Unpacking objects:  27% (5/18)   
Unpacking objects:  33% (6/18)   
Unpacking objects:  38% (7/18)   
Unpacking objects:  44% (8/18)   
Unpacking objects:  50% (9/18)   
Unpacking objects:  55% (10/18)   
Unpacking objects:  61% (11/18)   
Unpacking objects:  66% (12/18)   
Unpacking objects:  72% (13/18)   
Unpacking objects:  77% (14/18)   
Unpacking objects:  83% (15/18)   
Unpacking objects:  88% (16/18)   
Unpacking objects:  94% (17/18)   
Unpacking objects: 100% (18/18)   
Unpacking objects: 100% (18/18), done.
Checking connectivity... done.
einsteinium$ pwd
/tmp
einsteinium$ ls
#a#
A2CD1490-87C9-4C25-B83A-9BFD1C0CA4AC_IN
A2CD1490-87C9-4C25-B83A-9BFD1C0CA4AC_OUT
KSOutOfProcessFetcher.501.IW-ShwaXJwmHlkfc5wb2ZMHOaJY=
Makefile
UniMath
Workshop-2015-Boise  <---
a
com.apple.launchd.DCbIf7imEu
com.apple.launchd.TRXhLgKfCM
com.apple.launchd.cP0aa0C0Zf
emacs501
steam_chrome_shmem
einsteinium$ cd Workshop-2015-Boise/
einsteinium$ ls
Loci.m2		README.md	Splines
einsteinium$ ls -R
Loci.m2		README.md	Splines

./Splines:
foo.m2
einsteinium$ mkdir git-intro
einsteinium$ cd git-intro/
einsteinium$ ls
README.txt
einsteinium$ git add README.txt 
einsteinium$ git commit 
Waiting for Emacs...
[master 022d1ee] first commit of README.txt
 1 file changed, 1 insertion(+)
 create mode 100644 git-intro/README.txt
einsteinium$ git status
On branch master
Your branch is ahead of 'origin/master' by 1 commit.
  (use "git push" to publish your local commits)
nothing to commit, working directory clean
einsteinium$ git push
Username for 'https://github.com': DanGrayson
Password for 'https://DanGrayson@github.com': 
Counting objects: 4, done.
Delta compression using up to 8 threads.
Compressing objects:  50% (1/2)   
Compressing objects: 100% (2/2)   
Compressing objects: 100% (2/2), done.
Writing objects:  25% (1/4)   
Writing objects:  50% (2/4)   
Writing objects:  75% (3/4)   
Writing objects: 100% (4/4)   
Writing objects: 100% (4/4), 343 bytes | 0 bytes/s, done.
Total 4 (delta 1), reused 0 (delta 0)
To https://github.com/Macaulay2/Workshop-2015-Boise.git
   da55215..022d1ee  master -> master
einsteinium$ git status
On branch master
Your branch is up-to-date with 'origin/master'.
nothing to commit, working directory clean
einsteinium$ git pull
Already up-to-date.
einsteinium$ ls
README.txt
einsteinium$ git status
On branch master
Your branch is up-to-date with 'origin/master'.
Changes not staged for commit:
  (use "git add <file>..." to update what will be committed)
  (use "git checkout -- <file>..." to discard changes in working directory)

	modified:   README.txt

no changes added to commit (use "git add" and/or "git commit -a")
einsteinium$ git commit -a
Waiting for Emacs...
[master 64e67d4] added a line
 1 file changed, 2 insertions(+)
einsteinium$ git status
On branch master
Your branch is ahead of 'origin/master' by 1 commit.
  (use "git push" to publish your local commits)
nothing to commit, working directory clean
einsteinium$ git push
Username for 'https://github.com': DanGrayson
Password for 'https://DanGrayson@github.com': 
To https://github.com/Macaulay2/Workshop-2015-Boise.git
 ! [rejected]        master -> master (fetch first)
error: failed to push some refs to 'https://github.com/Macaulay2/Workshop-2015-Boise.git'
hint: Updates were rejected because the remote contains work that you do
hint: not have locally. This is usually caused by another repository pushing
hint: to the same ref. You may want to first integrate the remote changes
hint: (e.g., 'git pull ...') before pushing again.
hint: See the 'Note about fast-forwards' in 'git push --help' for details.
einsteinium$ git pull
remote: Counting objects: 8, done.        
remote: Compressing objects:  12% (1/8)           
remote: Compressing objects:  25% (2/8)           
remote: Compressing objects:  37% (3/8)           
remote: Compressing objects:  50% (4/8)           
remote: Compressing objects:  62% (5/8)           
remote: Compressing objects:  75% (6/8)           
remote: Compressing objects:  87% (7/8)           
remote: Compressing objects: 100% (8/8)           
remote: Compressing objects: 100% (8/8), done.        
remote: Total 8 (delta 2), reused 0 (delta 0), pack-reused 0        
Unpacking objects:  12% (1/8)   
Unpacking objects:  25% (2/8)   
Unpacking objects:  37% (3/8)   
Unpacking objects:  50% (4/8)   
Unpacking objects:  62% (5/8)   
Unpacking objects:  75% (6/8)   
Unpacking objects:  87% (7/8)   
Unpacking objects: 100% (8/8)   
Unpacking objects: 100% (8/8), done.
From https://github.com/Macaulay2/Workshop-2015-Boise
   022d1ee..6925999  master     -> origin/master
Waiting for Emacs...
Merge made by the 'recursive' strategy.
 Splines/Splines.m2 | 1042 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 1042 insertions(+)
 create mode 100644 Splines/Splines.m2
einsteinium$ git log -1
commit abea298d49b174db3a90b9230cd0d5688711143f
Merge: 64e67d4 6925999
Author: Daniel R. Grayson <dan@math.uiuc.edu>
Date:   Wed May 27 08:57:09 2015 -0600

    Merge branch 'master' of https://github.com/Macaulay2/Workshop-2015-Boise
einsteinium$ git log -2
commit abea298d49b174db3a90b9230cd0d5688711143f
Merge: 64e67d4 6925999
Author: Daniel R. Grayson <dan@math.uiuc.edu>
Date:   Wed May 27 08:57:09 2015 -0600

    Merge branch 'master' of https://github.com/Macaulay2/Workshop-2015-Boise

commit 64e67d4397010ad01e6d2d54c7f246c4642bde05
Author: Daniel R. Grayson <dan@math.uiuc.edu>
Date:   Wed May 27 08:56:00 2015 -0600

    added a line
einsteinium$ git log -3
commit abea298d49b174db3a90b9230cd0d5688711143f
Merge: 64e67d4 6925999
Author: Daniel R. Grayson <dan@math.uiuc.edu>
Date:   Wed May 27 08:57:09 2015 -0600

    Merge branch 'master' of https://github.com/Macaulay2/Workshop-2015-Boise

commit 64e67d4397010ad01e6d2d54c7f246c4642bde05
Author: Daniel R. Grayson <dan@math.uiuc.edu>
Date:   Wed May 27 08:56:00 2015 -0600

    added a line

commit 69259996859f55f4abda83c244d9a84c2b2009f7
Author: Gwyn Whieldon <whieldon@hood.edu>
Date:   Wed May 27 08:55:29 2015 -0600

    Added small change.
einsteinium$ git log -3 -u
commit abea298d49b174db3a90b9230cd0d5688711143f
Merge: 64e67d4 6925999
Author: Daniel R. Grayson <dan@math.uiuc.edu>
Date:   Wed May 27 08:57:09 2015 -0600

    Merge branch 'master' of https://github.com/Macaulay2/Workshop-2015-Boise

commit 64e67d4397010ad01e6d2d54c7f246c4642bde05
Author: Daniel R. Grayson <dan@math.uiuc.edu>
Date:   Wed May 27 08:56:00 2015 -0600

    added a line

diff --git a/git-intro/README.txt b/git-intro/README.txt
index fe8c8b2..aa20feb 100644
--- a/git-intro/README.txt
+++ b/git-intro/README.txt
@@ -1 +1,3 @@
 hi there!!!
+and hello!
+

commit 69259996859f55f4abda83c244d9a84c2b2009f7
Author: Gwyn Whieldon <whieldon@hood.edu>
Date:   Wed May 27 08:55:29 2015 -0600

    Added small change.

diff --git a/Splines/Splines.m2 b/Splines/Splines.m2
index a65ea33..56bef11 100644
--- a/Splines/Splines.m2
+++ b/Splines/Splines.m2
@@ -1,4 +1,5 @@
 --Package for computing topological boundary maps and piecewise continuous splines on polyhedral complexes  --
+--Making a small change
 
 newPackage("Splines",DebuggingMode => true)
 
einsteinium$ gitk
  C-c C-z
[1]+  Stopped                 gitk
einsteinium$ bg
[1]+ gitk &
einsteinium$ git push
Username for 'https://github.com': DanGrayson
Password for 'https://DanGrayson@github.com': 
remote: Invalid username or password.
fatal: Authentication failed for 'https://github.com/Macaulay2/Workshop-2015-Boise.git/'
einsteinium$ git push
Username for 'https://github.com': DanGrayson
Password for 'https://DanGrayson@github.com': 
To https://github.com/Macaulay2/Workshop-2015-Boise.git
 ! [rejected]        master -> master (fetch first)
error: failed to push some refs to 'https://github.com/Macaulay2/Workshop-2015-Boise.git'
hint: Updates were rejected because the remote contains work that you do
hint: not have locally. This is usually caused by another repository pushing
hint: to the same ref. You may want to first integrate the remote changes
hint: (e.g., 'git pull ...') before pushing again.
hint: See the 'Note about fast-forwards' in 'git push --help' for details.
einsteinium$ git pull
remote: Counting objects: 11, done.        
remote: Compressing objects:   9% (1/11)           
remote: Compressing objects:  18% (2/11)           
remote: Compressing objects:  27% (3/11)           
remote: Compressing objects:  36% (4/11)           
remote: Compressing objects:  45% (5/11)           
remote: Compressing objects:  54% (6/11)           
remote: Compressing objects:  63% (7/11)           
remote: Compressing objects:  72% (8/11)           
remote: Compressing objects:  81% (9/11)           
remote: Compressing objects:  90% (10/11)           
remote: Compressing objects: 100% (11/11)           
remote: Compressing objects: 100% (11/11), done.        
Unpacking objects:   9% (1/11)   
Unpacking objects:  18% (2/11)   
Unpacking objects:  27% (3/11)   
remote: Total 11 (delta 2), reused 0 (delta 0), pack-reused 0        
Unpacking objects:  36% (4/11)   
Unpacking objects:  45% (5/11)   
Unpacking objects:  54% (6/11)   
Unpacking objects:  63% (7/11)   
Unpacking objects:  72% (8/11)   
Unpacking objects:  81% (9/11)   
Unpacking objects:  90% (10/11)   
Unpacking objects: 100% (11/11)   
Unpacking objects: 100% (11/11), done.
From https://github.com/Macaulay2/Workshop-2015-Boise
   6925999..99fabe5  master     -> origin/master
Waiting for Emacs...
Merge made by the 'recursive' strategy.
 .../August2014ConformalBlocksCode/ConformalBlocks.m2         | 1486 +++++++++++++++++++++++++
 QuantumCohomology/August2014ConformalBlocksCode/LieTypes.m2  | 1446 ++++++++++++++++++++++++
 Splines/SplinesNew.m2                                        |    1 +
 Splines/{Splines.m2 => SplinesOriginal.m2}                   |    0
 4 files changed, 2933 insertions(+)
 create mode 100644 QuantumCohomology/August2014ConformalBlocksCode/ConformalBlocks.m2
 create mode 100644 QuantumCohomology/August2014ConformalBlocksCode/LieTypes.m2
 create mode 100644 Splines/SplinesNew.m2
 rename Splines/{Splines.m2 => SplinesOriginal.m2} (100%)
einsteinium$ git push
Username for 'https://github.com': DanGrayson
Password for 'https://DanGrayson@github.com': 
Counting objects: 8, done.
Delta compression using up to 8 threads.
Compressing objects:  16% (1/6)   
Compressing objects:  33% (2/6)   
Compressing objects:  50% (3/6)   
Compressing objects:  66% (4/6)   
Compressing objects:  83% (5/6)   
Compressing objects: 100% (6/6)   
Compressing objects: 100% (6/6), done.
Writing objects:  12% (1/8)   
Writing objects:  25% (2/8)   
Writing objects:  37% (3/8)   
Writing objects:  50% (4/8)   
Writing objects:  62% (5/8)   
Writing objects:  75% (6/8)   
Writing objects:  87% (7/8)   
Writing objects: 100% (8/8)   
Writing objects: 100% (8/8), 923 bytes | 0 bytes/s, done.
Total 8 (delta 3), reused 0 (delta 0)
To https://github.com/Macaulay2/Workshop-2015-Boise.git
   99fabe5..e6e5259  master -> master

