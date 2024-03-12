import os
import smtplib
from astropy.io import fits
from email.message import EmailMessage

DIFFERENCED_PATH = 'differenced/'
FLAGGED_PATH = 'flagged/'

def send_email_with_fits(filepath, recipient_email, sender_email, sender_password, smtp_server, smtp_port):
    # Directly define the subject and body of the email
    subject = "Anomaly Detected in FITS Image: {filename}".format(filename=os.path.basename(filepath))
    
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
            if key.startswith('XCENT'):
                xcentroids.append(header[key])
            elif key.startswith('YCENT'):
                ycentroids.append(header[key])

        # Construct coordinates string
        coordinates_str = "\n".join([f"Coordinates of anomaly {i+1} (x, y): ({x}, {y})"
                                     for i, (x, y) in enumerate(zip(xcentroids, ycentroids))])

        # Format the body with the extracted information
        body = """
        Name of the image file: {filename}
        Time image was taken: {time_image_taken}
        {coordinates}
        RA and Dec of the FITS image: (RA: {ra}, Dec: {dec})
        """.format(
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

    fn = "NEOS_SCI_2021317013021.fits"
                
    with fits.open(DIFFERENCED_PATH + fn, mode='update') as hdul:
        header = hdul[0].header
        header[f'XCENT_1'] = 1
        header[f'YCENT_1'] = 1
                
        # Save changes to header
        hdul.flush()  # Writes the updated header to the file
                
        # Now move the updated file to the flagged folder
        os.rename(DIFFERENCED_PATH + fn, FLAGGED_PATH + fn)
        print(f"Flagged {fn} for classification.")


   
    filepath = os.path.join(FLAGGED_PATH, fn)
    try:
        send_email_with_fits(
            filepath=filepath,
            recipient_email='amelie.barsoum@gmail.com', 
            sender_email='capstone.18.neossat@gmail.com',   
            sender_password='pbta yiit hati tvkj',    
            smtp_server='smtp.gmail.com',           
            smtp_port=465                             
        )
    except Exception as e:
        print(f"An error occurred while sending the email for {fn}: {e}")