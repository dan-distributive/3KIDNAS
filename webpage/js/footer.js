async function loadFooter(){
	const res=await fetch('/components/footer.html');
	const html=await res.text();
	document.getElementById('site-footer-placeholder').innerHTML=html;
	
  // Set current year
  document.querySelector('.year').textContent = new Date().getFullYear();
}

document.addEventListener('DOMContentLoaded', loadFooter);