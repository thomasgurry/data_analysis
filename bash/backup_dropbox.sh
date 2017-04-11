#!/bin/bash

# Script to backup Dropbox contents.

# NOTES:
#        Writes to ~/lockup/backups/dropbox/

today=`eval 'date +"%m-%d-%y"'`
backup_dir=~/lockup/tornam/backups/dropbox/dropbox_backup_$today

if [ ! -d "$backup_dir" ]; then
    mkdir $backup_dir
    chmod 755 $backup_dir
    echo "Backing up Dropbox to directory 'dropbox_backup_$today'.  Please wait."
    cp -r ~/Dropbox/* $backup_dir &
    wait
    echo "You successfully backed up your Dropbox.  Now fuck off."
else
    echo "Already backed up today."
fi

