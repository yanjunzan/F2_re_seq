#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> git tutorial>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# There are many basics, but here I skipped all of them
# first clone the folder I have uploaded
git clone https://github.com/yanjunzan/F2_re_seq.git
# you should see stuff in the folder

##There are many elegant way to avoid modify remote repository, i.e, there are something 
##called branch . but for now  the simple one is do not modify it until you are 100% sure. 
# Example on how to upload modified script
# do some thing on ./R/README_oc.R
git add ./R/README_oc.R ## after change a file, git add filename should be ran
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
