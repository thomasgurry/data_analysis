#!/bin/bash

# Script to backup Dropbox/Admin contents.

# NOTES:
#        Writes to ~/lockup/backups/admin/

today=`eval 'date +"%m-%d-%y"'`
backup_dir=~/lockup/tornam/backups/admin/admin_backup_$today

if [ ! -d "$backup_dir" ]; then
    mkdir $backup_dir
    chmod 755 $backup_dir
    echo "Backing up admin to directory 'admin_backup_$today'.  Please wait."
    cp -r ~/Dropbox/Admin/* $backup_dir &
    wait
    echo "You successfully backed up your admin.  Now fuck off."
else
    echo "Already backed up today."
fi

