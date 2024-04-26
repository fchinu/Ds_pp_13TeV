import os
import io
from googleapiclient.discovery import build
from google.oauth2.credentials import Credentials
from google.auth.transport.requests import Request
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.exceptions import RefreshError
from googleapiclient.http import MediaIoBaseUpload
from PIL import Image

# Define the scopes needed for Google Slides API and Google Drive API
SCOPES = ['https://www.googleapis.com/auth/drive', 'https://www.googleapis.com/auth/presentations']

def authenticate():
    """Authenticate with Google APIs."""
    creds = None
    if os.path.exists('token.json'):
        creds = Credentials.from_authorized_user_file('token.json')
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            try:
                creds.refresh(Request())
            except RefreshError:
                flow = InstalledAppFlow.from_client_secrets_file(
                    'credentials.json', SCOPES)
                creds = flow.run_local_server(port=0)
        else:
            flow = InstalledAppFlow.from_client_secrets_file(
                'credentials.json', SCOPES)
            creds = flow.run_local_server(port=0)
        with open('token.json', 'w') as token:
            token.write(creds.to_json())
    return creds

def create_slide(service, presentation_id, image_path, creds):
    """Create a slide with an image."""
    image = Image.open(image_path)
    image_bytes = io.BytesIO()
    image.save(image_bytes, format='PNG')
    image_content = image_bytes.getvalue()
    
    # Upload image to Drive
    drive_service = build('drive', 'v3', credentials=creds)
    file_metadata = {
        'name': 'temp_image1.png',
        'parents': ['root'],
        'mimeType': 'image/png',
    }
    media = MediaIoBaseUpload(io.BytesIO(image_content), mimetype='image/png')
    drive_response = drive_service.files().create(
        body=file_metadata,
        media_body=media).execute()
    image_drive_id = drive_response.get('id')

    # Make the file publicly accessible
    drive_service.permissions().create(
        fileId=image_drive_id,
        body={'role': 'reader', 'type': 'anyone'}).execute()

    # Get the web link for the image
    rsp = drive_service.files().list(q=f"name='{file_metadata['name']}'").execute()['files'][0]
    print(rsp)
    print(f"Found file{rsp['name']} with id {rsp['id']}")
    image_web_link = f"{drive_service.files().get_media(fileId=rsp['id']).uri}&access_token={creds.token}"
    #image_web_link = drive_service.files().get_media(fileId=file_metadata['name']).uri + f"&access_token={creds.token}"
    print(image_web_link)

    # Add image to slide
    requests = [
        {
            'createImage': {
                'url': image_web_link,
                'elementProperties': {
                    'pageObjectId': 'g1',  # ID of the slide
                    'size': {
                        'height': {
                            'magnitude': image.height * 9525,
                            'unit': 'EMU'
                        },
                        'width': {
                            'magnitude': image.width * 9525,
                            'unit': 'EMU'
                        }
                    },
                    'transform': {
                        'scaleX': 1,
                        'scaleY': 1,
                        'translateX': 100000,
                        'translateY': 100000,
                        'unit': 'EMU'
                    }
                }
            }
        }
    ]

    # Execute the requests
    service.presentations().batchUpdate(
        presentationId=presentation_id, body={'requests': requests}).execute()

    # Delete the image from Drive
    drive_service.files().delete(fileId=image_drive_id).execute()


def main():
    # Authenticate with Google APIs
    creds = authenticate()

    # Build the service
    service = build('slides', 'v1', credentials=creds)

    # Create a new presentation
    presentation = service.presentations().create(
        body={'title': 'My Presentation'}).execute()
    presentation_id = presentation.get('presentationId')

    # Path to the image file
    image_path = '/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/figure.png'

    # Create a slide with the image
    create_slide(service, presentation_id, image_path, creds)

    print(f'Presentation created: https://docs.google.com/presentation/d/{presentation_id}')

if __name__ == '__main__':
    main()
