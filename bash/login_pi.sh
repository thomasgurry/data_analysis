#!/usr/bin/expect

set login "pi"
set addr "67.186.134.9"
set pw "gasberry"

spawn ssh $login@$addr
expect "$login@$addr\'s password:"
send "$pw\r"
interact
