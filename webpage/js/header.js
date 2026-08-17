async function loadHeader(){
  const res = await fetch('/components/header.html');
  const html = await res.text();
  document.getElementById('site-header-placeholder').innerHTML = html;
  
  const path = window.location.pathname.split('/').pop() || 'index.html';
  const keyMap = {
    'platform.html':'platform',
    'suitability.html':'suitability',
    'security-job-data.html':'jobs',
    'security-worker.html':'workers',
    'workers.html':'workers',
    'worker-browser.html':'workers',
    'worker-docker.html':'workers',
    'worker-linux.html':'workers',
    'worker-windows.html':'workers',
    'compute-economics.html':'economics',
    'about.html':'about',
    'privacy-policy.html':'about',
    'index.html':'platform'
  };
  const key = keyMap[path] || '';
  document.querySelectorAll('.nav-item').forEach(n=>{
    if(n.dataset.key === key){ n.classList.add('active'); }
  });
  
  const portal = document.getElementById('portalBtn');
  if(portal){
    portal.addEventListener('click', ()=>{
      portal.classList.add('clicked');
      setTimeout(()=> portal.classList.remove('clicked'), 300);
      window.location.href = 'portal.html';
    });
  }
  
  // *** HAMBURGER MENU CODE - ADD THIS HERE ***
  const hamburger = document.getElementById('hamburger');
  const navList = document.getElementById('navList');
  
  console.log('Hamburger:', hamburger);
  console.log('NavList:', navList);
  
  if (hamburger && navList) {
    hamburger.addEventListener('click', function(e) {
      e.preventDefault();
      console.log('Hamburger clicked!');
      hamburger.classList.toggle('active');
      navList.classList.toggle('active');
      console.log('NavList classes:', navList.classList);
      console.log('Hamburger classes:', hamburger.classList);
    });
    
    // Close menu when clicking a nav item
    const navItems = document.querySelectorAll('.nav-item');
    navItems.forEach(item => {
      item.addEventListener('click', () => {
        hamburger.classList.remove('active');
        navList.classList.remove('active');
      });
    });
  } else {
    console.error('Hamburger or NavList not found!');
  }
  // *** END HAMBURGER MENU CODE ***
}

document.addEventListener('DOMContentLoaded', loadHeader);

document.addEventListener("DOMContentLoaded", () => {
  const dropdown = document.querySelector(".dropdown");
  if (dropdown) {
    const aboutLink = dropdown.querySelector(".nav-item");
    if (aboutLink) {
      // Toggle open on click (mobile-friendly)
      aboutLink.addEventListener("click", (e) => {
        // On mobile, prevent navigating immediately
        if (window.innerWidth < 900) {
          e.preventDefault();
          dropdown.classList.toggle("open");
        }
      });
    }
  }
});

function copyToClipboard(text) {
  // Modern approach using Clipboard API
  if (navigator.clipboard && navigator.clipboard.writeText) {
    navigator.clipboard.writeText(text).then(() => {
      console.log('Email copied to clipboard!');
      // Optional: Show a temporary success message
      showCopyFeedback();
    }).catch(err => {
      console.error('Failed to copy:', err);
    });
  } else {
    // Fallback for older browsers
    const textArea = document.createElement('textarea');
    textArea.value = text;
    textArea.style.position = 'fixed';
    textArea.style.left = '-999999px';
    document.body.appendChild(textArea);
    textArea.select();
    try {
      document.execCommand('copy');
      console.log('Email copied to clipboard!');
      showCopyFeedback();
    } catch (err) {
      console.error('Failed to copy:', err);
    }
    document.body.removeChild(textArea);
  }
}

// Optional: Visual feedback function
function showCopyFeedback() {
  const button = event.target.closest('button');
  const originalText = button.querySelector('.btn-text').textContent;
  button.querySelector('.btn-text').textContent = 'COPIED!';
  
  setTimeout(() => {
    button.querySelector('.btn-text').textContent = originalText;
  }, 2000);
}
