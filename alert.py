#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 23:42:10 2022

@author: luischacon
"""


def alert(subject,body,to):
    import smtplib
    from email.message import EmailMessage
    
    user = 'lch285@g.uky.edu'
    password = 'vhjacehswgdbnsck'
    msg = EmailMessage()
    msg.set_content(body)
    msg['subject'] = subject
    msg['to'] = to
    msg['from'] = user
    
    
    server = smtplib.SMTP('smtp.gmail.com',587)
    server.starttls()
    
    server.login(user, password)
    server.send_message(msg)
    server.quit()


if __name__ == '__main__':
    
    email = 'lch285@g.uky.edu'
    phone = '8594901117@txt.att.net' # works for ATT
    
    method = 'text'
    # method = 'email'
    
    if method == 'text':
        to = phone
    else:
        to = email
   
    alert('Test', 'This is a test', to)
