#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> git tutorial>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# There are many basics, but here I skipped all of them
# first clone the folder I have uploaded
git clone https://github.com/yanjunzan/F2_re_seq.git
# you should see folders and code downloaded (you will find the updated code here)

##There are many elegant way to avoid modify remote repository, i.e, there are something 
##called branch . but for now  the simple one is do not modify it until you are 100% sure. 
# Example on how to upload modified script
# Imagine that some thing have been changed for ./R/README_oc.R
git add ./R/README_oc.R ##  git add filename , tell git I have changed this 
# git diff will show what has been changed.
git commit -m "I changed this" # this will confirm the change and  need to write some note
# push this change to remote
git push

## if it does not work you mihgt need to set this up. I have invite you to this project
## please confirm it first
git remote add origin https://github.com/yanjunzan/F2_re_seq.git

## for this git push to run, a few configuration needs to set
# id stuff
 git config --global user.name "orjancarlborg" # this is your id, you sent to me ages ago
 git config --global user.email yanjunzan@gmail.com # change to your email
# default editor . change vim to the text editor you comfortable with. e.g. emacs ..
 git config --global core.editor vim
# double check he set up
 git config --list ### it list the file located in /etc/gitconfig or  ~/.gitconfig

 
 ########>>>>>>>>>>>>>>>>>>>>>>>>>>>> 171215. git handling conflit
 #######<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<branching<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 # you do not want to always overwrite the files in remote server, so need a backup files called branch, 
 #everyone is working on their own branch and later merged with master(server file)
 
 # now create a branch
 $ git branch yanjun
 # switch to the new branch
 $ git checkout yanjun
 $ git branch # this will shown wich branch you are currently in 
 # do some changes, and commit it by
 $ git commit -m "changed sth" #everything changed here, in only visible to me in thsi branch. 
 
 ## After I am done with my changes, I want to merge it with master and push back to server
 ######### >>>>>>>>>>>>>>>>>>>>>>>>>>>merge branching<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 $ git checkout master
 $ git merge yanjun 
 $ git push # now changes in branch yanjun is already merged with master, this command will overwrite
 ## remote server. and conflicts will be marked and reminded.
 ## delete my branch, you do not have to do this
 $ git branch -D yanjun
 
 #### there are more on how to push branch to server, we will learn it later
 