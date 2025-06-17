import smtplib
import csv
import re
from email.message import EmailMessage
from email.utils import make_msgid
from pathlib import Path

# --- Configuration ---
smtp_server = "" #### smtp server
smtp_port = 25 
sender_email = "fdumetz@som.umaryland.edu" ### email address

# --- Helper: escape CSS curly braces but keep {name}
def escape_curly_braces_except_placeholders(text, keep_keys):
    escaped = text.replace("{", "{{").replace("}", "}}")
    for key in keep_keys:
        escaped = escaped.replace("{{" + key + "}}", "{" + key + "}")
    return escaped

# --- Load contacts ---
contacts = []
csv_path = "contacts.csv"   ### see contact.csv example

with open(csv_path, newline='', encoding="utf-8-sig") as csvfile:
    reader = csv.DictReader(csvfile, delimiter=',')
    reader.fieldnames = [h.strip().lower() for h in reader.fieldnames]
    print("Detected headers:", reader.fieldnames)

    for row in reader:
        row = {k.strip().lower(): v for k, v in row.items()}
        name = row.get("name")
        email = row.get("email")
        if name and email:
            contacts.append((name.strip(), email.strip()))

print(f"✅ Loaded {len(contacts)} contacts.")

# --- Load and escape HTML template ---
template_path = "email_template.html" ### make an html file on word
with open(template_path, "r", encoding="utf-8", errors="replace") as f:
    raw_html = f.read()

# Escape CSS curly braces but keep {name}
html_template = escape_curly_braces_except_placeholders(raw_html, keep_keys=["name"])

# --- Replace local image paths with CID links ---
image_folder = Path("email_template.fld")
image_files = re.findall(r'src="email_template\.fld/([^"]+)"', html_template)
cid_map = {img: make_msgid(domain="igs.umaryland.edu")[1:-1] for img in image_files}

for img, cid in cid_map.items():
    html_template = html_template.replace(f'email_template.fld/{img}', f'cid:{cid}')

# --- Start SMTP connection ---
server = smtplib.SMTP(smtp_server, smtp_port)
server.set_debuglevel(1)

# --- Send personalized emails ---
for name, email in contacts:
    print(f"\nPreparing message for: {name} <{email}>")

    msg = EmailMessage()
    msg["From"] = sender_email
    msg["To"] = email
    msg["Subject"] = "Reminder to renew your ASTMH Membership"
    msg["Reply-To"] = sender_email

    personalized_html = html_template.format(name=name)
    plain_fallback = f"Hi {name},\nThis is an HTML email. Please view in a compatible email client."

    msg.set_content(plain_fallback)
    msg.add_alternative(personalized_html, subtype='html')

    # Embed inline images
    for img_name, cid in cid_map.items():
        img_path = image_folder / img_name
        if img_path.exists():
            with open(img_path, 'rb') as img:
                msg.get_payload()[1].add_related(
                    img.read(),
                    maintype='image',
                    subtype=img_path.suffix.lstrip('.'),
                    cid=f"<{cid}>"
                )
        else:
            print(f"⚠️ Image not found: {img_path}")

    server.send_message(msg)
    print(f"✅ Email sent to {email}")

server.quit()
print("🎉 All emails sent.")
