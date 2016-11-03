import logging.handlers
import smtplib
from threading import Thread


def smtp_at_your_own_leasure(mailhost, port, username, password, fromaddr, toaddrs, msg):
    smtp = smtplib.SMTP(mailhost, port)
    if username:
        smtp.ehlo() # for tls add this line
        smtp.starttls() # for tls add this line
        smtp.ehlo() # for tls add this line
        smtp.login(username, password)
    smtp.sendmail(fromaddr, toaddrs, msg)
    smtp.quit()


class ThreadedTlsSMTPHandler(logging.handlers.SMTPHandler):
    def emit(self, record):
        try:
            import string # for tls add this line
            try:
                from email.utils import formatdate
            except ImportError:
                formatdate = self.date_time
            port = self.mailport
            if not port:
                port = smtplib.SMTP_PORT
            msg = self.format(record)
            msg = "From: %s\r\nTo: %s\r\nSubject: %s\r\nDate: %s\r\n\r\n%s" % (
                            self.fromaddr,
                            ','.join(self.toaddrs),
                            self.getSubject(record),
                            formatdate(), msg)
            thread = Thread(target=smtp_at_your_own_leasure, args=(self.mailhost, port, self.username, self.password, self.fromaddr, self.toaddrs, msg))
            thread.start()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)
