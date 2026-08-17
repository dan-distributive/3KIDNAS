document.addEventListener('DOMContentLoaded', () => {
  const content = document.getElementById('docContent');
  const tocList = document.getElementById('tocList');
  const headings = content.querySelectorAll('h2, h3');

  let currentH2Item = null;
  let currentH2Link = null;
  const headingMap = new Map();

  headings.forEach((heading, index) => {
    if (!heading.id) {
      heading.id = heading.textContent
        .toLowerCase()
        .replace(/[^a-z0-9]+/g, '-')
        .replace(/(^-|-$)/g, '') + '-' + index;
    }

    /* ---------- H2 ---------- */
    if (heading.tagName === 'H2') {
      const li = document.createElement('li');
      li.classList.add('toc-item');

      const a = document.createElement('a');
      a.href = `#${heading.id}`;
      a.textContent = heading.textContent;
      a.classList.add('toc-h2');

      const sublist = document.createElement('ul');
      sublist.classList.add('toc-sublist');

      a.addEventListener('click', e => {
        e.preventDefault();

        /* Close all other sections */
        document.querySelectorAll('.toc-item.open').forEach(item => {
          if (item !== li) item.classList.remove('open');
        });

        /* Toggle this one */
        if (!li.classList.contains('no-children')) {
          li.classList.toggle('open');
        }

        document.getElementById(heading.id).scrollIntoView({ behavior: 'smooth' });
      });

      li.appendChild(a);
      li.appendChild(sublist);
      tocList.appendChild(li);

      currentH2Item = sublist;
      currentH2Link = li;

      headingMap.set(heading.id, a);
    }

    /* ---------- H3 ---------- */
    if (heading.tagName === 'H3' && currentH2Item) {
      const subLi = document.createElement('li');
      const subA = document.createElement('a');

      subA.href = `#${heading.id}`;
      subA.textContent = heading.textContent;
      subA.classList.add('toc-h3');

      subLi.appendChild(subA);
      currentH2Item.appendChild(subLi);

      headingMap.set(heading.id, subA);
    }
  });

  /* ---------- Remove arrows from H2 with no H3 ---------- */
  document.querySelectorAll('.toc-item').forEach(item => {
    const sublist = item.querySelector('.toc-sublist');
    if (!sublist || sublist.children.length === 0) {
      item.classList.add('no-children');
    }
  });

  /* ---------- SCROLL SPY (highlight ONLY) ---------- */

  const observer = new IntersectionObserver(
    entries => {
      entries.forEach(entry => {
        if (entry.isIntersecting) {
          headingMap.forEach(link => link.classList.remove('active'));

          const activeLink = headingMap.get(entry.target.id);
          if (activeLink) {
            activeLink.classList.add('active');

            const parentItem = activeLink.closest('.toc-item');
            if (parentItem) {
              parentItem.querySelector('.toc-h2')?.classList.add('active');
            }
          }
        }
      });
    },
    {
      rootMargin: '-30% 0px -60% 0px',
      threshold: 0
    }
  );

  headings.forEach(h => observer.observe(h));
});