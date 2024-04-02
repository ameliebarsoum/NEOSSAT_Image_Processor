import os
import smtplib
from astropy.io import fits
from email.message import EmailMessage
import shutil

DIFFERENCED_PATH = '../ObjectDetection/differenced/'
FLAGGED_PATH = '../ObjectDetection/flagged_original_fits/'

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
        ras = []
        decs = []

        # Extract all XCENTER and YCENTER values
        for key in header.keys():
            if key.startswith('XCENT'):
                xcentroids.append(header[key])
            elif key.startswith('YCENT'):
                ycentroids.append(header[key])
            elif key.startswith('DETRA'):
                ras.append(header[key])
            elif key.startswith('DETDEC'):
                decs.append(header[key])
        
        print(xcentroids, ycentroids, ras, decs)
        # Construct coordinates string for the email body:
                # Format: Coordinates of source with id: id = (X: x, Y: y) -> (RA: ra, Dec: dec)
        coordinates_str = ""
        i = 1
        for x, y, ra, dec in zip(xcentroids, ycentroids, ras, decs):
            coordinates_str += (f"Coordinates of anomaly {i} = (X: {x}, Y: {y}) -> (RA: {ra}, Dec: {dec})\n")
            i += 1

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
    
    for fn in os.listdir(FLAGGED_PATH):

        filepath = os.path.join(FLAGGED_PATH, fn)
        try:
            send_email_with_fits(
                filepath=filepath,
                recipient_email='amelie.barsoum@gmail.com', 
                sender_email='capstone.18.neossat@gmail.com',   
                sender_password='melj eluf vndc gpwt',    
                smtp_server='smtp.gmail.com',           
                smtp_port=465                             
            )
            break
        except Exception as e:
            print(f"An error occurred while sending the email for {fn}: {e}")
            continue