import requests 
import os
BASE_URL = 'http://localhost:8002'

def request_access(name, surname, email, github_username, url=BASE_URL):
    headers = {
        "Content-Type": "application/json"
    }
    data = {
        "name": name,
        "surname": surname,
        "email": email,
        "github_username": github_username
    }
    response = requests.post(f'{url}/access/request', json=data, headers=headers)
    print(f"Status Code: {response.status_code}")
    print(f"Response: {response.text}")
    return response

def request_access_interactive(url=BASE_URL):
    name = input("Enter your name: ")
    surname = input("Enter your surname: ")
    email = input("Enter your email: ")
    github_username = input("Enter your GitHub username: ")
    
    return request_access(name, surname, email, github_username, url)

def get_server_token(url=BASE_URL):
    headers = {
        "Content-Type": "application/json"
    }
    
    # Check if GitHub token is available
    token = os.getenv('GITHUB_TOKEN')
    if not token:
        print("No token set. Please set your GitHub token as env variable with name GITHUB_TOKEN")
        return None
    
    data = {
        "token": token,
        "issue_token": True
    }
    
    try:
        response = requests.post(f'{url}/auth/github', json=data, headers=headers)
        
        response.raise_for_status()

        return response.json()["token"] 
        
    except requests.exceptions.ConnectionError:
        print(f"Error: Could not connect to server at {url}")
        return None
    except requests.exceptions.Timeout:
        print("Error: Request timed out")
        return None
    except requests.exceptions.HTTPError as e:
        print(f"Error: HTTP error occurred: {e}")
        print(f"Status code: {response.status_code}")
        if response.text:
            print(f"Response: {response.text}")
        return None
    except requests.exceptions.RequestException as e:
        print(f"Error: An error occurred while making the request: {e}")
        return None
    except Exception as e:
        print(f"Unexpected error: {e}")
        return None

def create_authenticated_session():
    session = requests.Session()
    token = get_server_token()
    if not token:
        print("Authentication Failed")
        return None
    session.headers.update({"Authorization": f"Bearer {token}"})
    print("Authentication Successfull")
    return session
    

