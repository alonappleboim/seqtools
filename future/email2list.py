#! /bin/python3

from datetime import datetime
import time
import csv
import argparse
import smtplib
import sys

SLEEP_INTERVAL = 300  # seconds
DEF_HOST = 'smtp.gmail.com'
DEF_PORT = 587
DEF_FROM = 'lab.chores@gmail.com'
DEF_USER = 'lab.chores'
DEF_PWD = 'labchores1234'


def send_email(transaction):
    t = transaction
    sess = smtplib.SMTP(t['host'], t['port'])
    sess.ehlo()
    sess.starttls()
    sess.login(t['user'], t['password'])
    text = 'From: %s\nTo: %s\nSubject: %s\n\n%s' % (t['fromaddr'], ", ".join(t['toaddrs']), t['subject'], t['msg'])
    sess.sendmail(t['fromaddr'], t['toaddrs'], text)
    sess.quit()


def parse_msg_list(file, args):
    commq = []
    for line, comm in enumerate(csv.DictReader(open(file, 'rU'), delimiter='\t')):
        to = comm['To'].split(',')
        try:
            date = comm['Date'].split('/')
            date = dict(day=int(date[0]), month=int(date[1]),
                        year=int(date[2]) if int(date[2]) > 2000 else int(date[2])+2000)
            date = datetime(**date)
        except TypeError:
            raise TypeError('date should be in "dd/mm/yyyy" format')
        try:
            msg = comm['Message']
        except KeyError:
            if args.default_message is None:
                raise ValueError('Default message not provided, and missing from line %i in message file' % line)
            msg = args.default_message
        try:
            subject= comm['Subject']
        except KeyError:
            if args.default_subject is None:
                raise ValueError('Default subject not provided, and missing from line %i in message file' % line)
            subject = args.default_subject
        commq.append(dict(toaddrs=to, date=date, msg=msg, subject=subject))
    return sorted(commq, key=lambda x: x['date'])


def main(commq, args):
    for comm in commq:
        while True:
            if not comm['date'].date() >= datetime.now().date(): time.sleep(SLEEP_INTERVAL)
            else:
                t = dict(args.__dict__)
                t.update(comm)
                send_email(t)
                if args.verbose:
                    sys.stderr.write('sent email to %s at %s.\n' % (','.join(comm['toaddrs']), str(datetime.now())))
                break


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('date_address_tsv_file', type=str,
                   help=('A tab-separated file with fields: To - comma separated email list; Date - a dd/mm/yyyy format on '
                         'which email will be sent; Msg - optional message field, otherwise --default_message is used, '
                         'Subject - optional subject field, otherwise --default_subject is used'))
    p.add_argument('--host', '-H', type=str, default=DEF_HOST, help='sender smtp host')
    p.add_argument('--default_message', '-m', type=str, default=None,
                   help='message sent in all emails, unless specifiecd in csv file')
    p.add_argument('--default_subject', '-s', type=str, default=None,
                   help='subject used in all emails, unless specifiecd in csv file')
    p.add_argument('--port', '-p', type=int, default=DEF_PORT, help='smtp host login port')
    p.add_argument('--user', '-u', type=str, default=DEF_USER, help='username used to login into host server')
    p.add_argument('--password', '-w', type=str, default=DEF_PWD, help="user's password")
    p.add_argument('--fromaddr', '-f', type=str, default=DEF_FROM, help='email address from which to send messages')
    p.add_argument('--nverbose', '-nv', dest='verbose', action='store_false',
                   help='do not output processing info to stderr')
    args = p.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    main(parse_msg_list(args.date_address_tsv_file, args), args)