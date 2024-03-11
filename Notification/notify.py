import os
import smtplib
from astropy.io import fits
from email.message import EmailMessage
import configparser

FLAGGED_PATH = 'flagged/'

def send_email_with_fits(filepath, recipient_email, sender_email, sender_password, smtp_server, smtp_port):
    # Load the email template from the config file
    config = configparser.ConfigParser()
    config.read('email_template.ini')
    subject_template = config['email']['subject']
    body_template = config['email']['body']
    
    # Open the FITS file and extract header information
    with fits.open(filepath) as hdul:
        header = hdul[0].header
        time_image_taken = header['DATE-OBS']
        ra = header['RA']
        dec = header['DEC']

        # Initialize variables to hold coordinates
        xcentroids = []
        ycentroids = []

        # Extract all XCENTER and YCENTER values
        for key in header.keys():
            if key.startswith('XCENTER'):
                xcentroids.append(header[key])
            elif key.startswith('YCENTER'):
                ycentroids.append(header[key])

        # Construct coordinates string
        coordinates_str = "\n".join([f"Coordinates of anomaly {i+1} (x, y): ({x}, {y})"
                                     for i, (x, y) in enumerate(zip(xcentroids, ycentroids))])

    # Format the subject and body with the extracted information
    subject = subject_template.format(filename=os.path.basename(filepath))
    body = body_template.format(
        filename=os.path.basename(filepath),
        time_image_taken=time_image_taken,
        coordinates=coordinates_str,
        ra=ra,
        dec=dec
    )

    # Create the email message
    msg = EmailMessage()
    msg['Subject'] = subject
    msg['From'] = sender_email
    msg['To'] = recipient_email
    msg.set_content(body)

    # Attach the FITS file
    with open(filepath, 'rb') as f:
        file_data = f.read()
        file_name = os.path.basename(filepath)
        msg.add_attachment(file_data, maintype='application', subtype='fits', filename=file_name)

    # Send the email via SMTP
    with smtplib.SMTP_SSL(smtp_server, smtp_port) as server:
        server.login(sender_email, sender_password)
        server.send_message(msg)

    print("Email sent successfully!")


if __name__ == "__main__":
    # Loop through all FITS files in the FLAGGED_PATH directory
    for filename in os.listdir(FLAGGED_PATH):
        if filename.lower().endswith('.fits'):
            filepath = os.path.join(FLAGGED_PATH, filename)
            try:
                # Send an email for each FITS file
                send_email_with_fits(
                    filepath=filepath,
                    recipient_email='amelie.barsoum@gmail.com', 
                    sender_email='capstone.18.neossat@gmail.com',   
                    sender_password='Capstone18!',    
                    smtp_server='smtp.gmail.com',           
                    smtp_port=465                             
                )
            except Exception as e:
                print(f"An error occurred while sending the email for {filename}: {e}")