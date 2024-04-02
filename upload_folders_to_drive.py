from google_auth_oauthlib.flow import InstalledAppFlow
from googleapiclient.discovery import build
from googleapiclient.http import MediaFileUpload
from google.auth.transport.requests import Request
import os
import pickle
import os
import sys

NEOSSAT_folder_id = '1vIRlgqbucJPp8XnrYxRG6gUdvr8fzvMH'
# Scopes define the level of access you need
SCOPES = ['https://www.googleapis.com/auth/drive.file']

def login_to_google_drive():
    creds = None
    # The file token.pickle stores the user's access and refresh tokens, and is
    # created automatically when the authorization flow completes for the first
    # time.
    if os.path.exists('token.pickle'):
        with open('token.pickle', 'rb') as token:
            creds = pickle.load(token)
    # If there are no (valid) credentials available, let the user log in.
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file(
                'client_secret_392747387582-9r0v5b3c2cqgb28ukfql0o18i7g5c1ko.apps.googleusercontent.com.json', SCOPES)
            creds = flow.run_local_server(port=0)
        # Save the credentials for the next run
        with open('token.pickle', 'wb') as token:
            pickle.dump(creds, token)

    service = build('drive', 'v3', credentials=creds)
    return service

def upload_file(service, filename, path, mimetype='application/octet-stream', parent_id=None):
    """
    Uploads a file to Google Drive into the specified folder.
    
    Args:
    - service: Authenticated Google Drive service instance.
    - filename: Name of the file to create in Google Drive.
    - path: Path to the file on local system.
    - mimetype: MIME type of the file.
    - parent_id: The ID of the parent folder to upload the file into.
    """
    file_metadata = {'name': filename}
    if parent_id:
        file_metadata['parents'] = [parent_id]
    media = MediaFileUpload(path, mimetype=mimetype)
    file = service.files().create(body=file_metadata, media_body=media, fields='id').execute()
    # print(f'Uploaded file ID: {file.get("id")}')


def create_drive_folder(service, name, parent_id=None):
    """
    Create a folder on Google Drive.
    
    Args:
    - service: Authenticated Google Drive service instance.
    - name: The name of the folder.
    - parent_id: The ID of the parent folder under which to create the new folder.
    
    Returns:
    The ID of the newly created folder.
    """
    file_metadata = {
        'name': name,
        'mimeType': 'application/vnd.google-apps.folder'
    }
    if parent_id:
        file_metadata['parents'] = [parent_id]
    folder = service.files().create(body=file_metadata, fields='id').execute()
    return folder.get('id')


def upload_folder(service, folder_path, drive_folder_name, parent_id=None):
    """
    Uploads an entire folder to Google Drive under a specified parent folder.
    
    Args:
    - service: Authenticated Google Drive service instance.
    - folder_path: Local path to the folder to be uploaded.
    - drive_folder_name: The name of the folder to create on Google Drive.
    - parent_id: The ID of the parent folder under which to create the new folder (existing folder ID).
    """
    # This will now create a new folder under the existing folder specified by parent_id
    drive_folder_id = create_drive_folder(service, drive_folder_name, parent_id)
    for filename in os.listdir(folder_path):
        filepath = os.path.join(folder_path, filename)
        if os.path.isfile(filepath):
            upload_file(service, filename, filepath, parent_id=drive_folder_id)
        elif os.path.isdir(filepath):
            upload_folder(service, filepath, filename, parent_id=drive_folder_id)
        else:
            print(f'Error: {filename} is neither a file nor a folder. Skipped.')
    print(f'All files uploaded to folder "{drive_folder_name}" within the existing Drive folder.')


def delete_drive_folder(service, folder_id):
    """
    Deletes a folder or file on Google Drive by ID.

    Args:
    - service: Authenticated Google Drive service instance.
    - folder_id: The ID of the folder or file to delete.
    """
    service.files().delete(fileId=folder_id).execute()
    print(f'Folder/File with ID {folder_id} has been deleted.')

def find_drive_folder(service, name, parent_id):
    """
    Finds a folder on Google Drive by name under a specified parent folder.

    Args:
    - service: Authenticated Google Drive service instance.
    - name: The name of the folder to find.
    - parent_id: The ID of the parent folder to search within.

    Returns:
    The ID of the found folder, or None if not found.
    """
    query = f"name='{name}' and '{parent_id}' in parents and mimeType='application/vnd.google-apps.folder' and trashed=false"
    response = service.files().list(q=query, spaces='drive', fields='files(id, name)').execute()
    files = response.get('files', [])
    if len(files) > 0:
        return files[0].get('id')
    return None

def create_or_recreate_drive_folder(service, name, parent_id):
    """
    Creates or recreates a folder on Google Drive. If the folder already exists,
    it is deleted and a new one is created.

    Args:
    - service: Authenticated Google Drive service instance.
    - name: The name of the folder.
    - parent_id: The ID of the parent folder under which to create the new folder.

    Returns:
    The ID of the newly created folder.
    """
    existing_folder_id = find_drive_folder(service, name, parent_id)
    if existing_folder_id:
        delete_drive_folder(service, existing_folder_id)
    return create_drive_folder(service, name, parent_id)


if __name__ == '__main__':

    args = sys.argv
    parent_folder = args[1]
    parent_folder = parent_folder.split('/')[-1]    # Extract the last folder name
    
    top_folder_id = NEOSSAT_folder_id    

    service = login_to_google_drive()
    parent_id = None

    flagged_difference_path = "ObjectDetection/flagged"
    if os.path.exists(flagged_difference_path) and len(os.listdir(flagged_difference_path)) > 0:
        if parent_id == None:
            parent_id = create_or_recreate_drive_folder(service, parent_folder, top_folder_id)
        upload_folder(service, flagged_difference_path, "Flagged_differenced", parent_id=parent_id)

    flagged_original_path = "ObjectDetection/flagged_original"
    if os.path.exists(flagged_original_path) and len(os.listdir(flagged_original_path)) > 0:
        if parent_id == None:
            parent_id = create_or_recreate_drive_folder(service, parent_folder, top_folder_id)
        upload_folder(service, flagged_original_path, "Flagged_original", parent_id=parent_id)

    cosmic_ray_hits_path = "ObjectDetection/cosmic_ray_hits"
    if os.path.exists(cosmic_ray_hits_path) and len(os.listdir(cosmic_ray_hits_path)) > 0:
        if parent_id == None:
            parent_id = create_or_recreate_drive_folder(service, parent_folder, top_folder_id)
        upload_folder(service, "ObjectDetection/cosmic_ray_hits", "Cosmic_ray_hits", parent_id=parent_id)
