/**
 *  @file       dcp-client.js
 *              Interface/loader for web browser consumers of the dcp-client bundle.
 *
 *              Once this script has been loaded in a vanilla-web environment, the dcp and dcpConfig
 *              symbols should be available to any scripts below this one in the DOM tree. This
 *              script is also aware of the following attributes on its own SCRIPT tag:
 *              - onload:    code to run when completely loaded; dcpConfig will be defined.
 *              - onready:   same as onload, but dcp is also defined.
 *              - scheduler: set the URL for the scheduler, also used to locate dcpConfig
 *              - bundle:    set the URL for the bundle to load, default is the one in this directory
 *  @author     Wes Garland, wes@kingsds.network
 *  @date       Aug 2019
 */
'use strict';
(function namespaceIIFE() {
  console.log(`%c
   _____ _____ ___________   _
  /  ___|_   _|  _  | ___ \\ | |
  \\ \`--.  | | | | | | |_/ / | |
   \`--. \\ | | | | | |  __/  |_|
  /\\__/ / | | \\ \\_/ / |      _
  \\____/  \\_/  \\___/\\_|     |_|

%c
The console is a browser feature intended for developers. If somebody told you to paste something here it may be a scam and your information could be stolen. Help us keep security in mind and keep your data safe.
  ~ DCP Team

https://distributive.network/`, "font-weight: bold; font-size: 1.2em; color: #00a473;", "font-size: 1.2em;");

  var _dcpConfig = typeof dcpConfig === 'object' ? dcpConfig : undefined;
  {
    let allScripts = document.getElementsByTagName('SCRIPT');
    let thisScript = allScripts[allScripts.length - 1];
    let thisScriptURL = new URL(thisScript.src)
    let schedulerURL;
    let dcpConfigHref = thisScript.getAttribute('dcpConfig') /* camel case is wrong+deprecated */
                     || thisScript.getAttribute('dcp-config');
    let configScript;

    /* kebab to camel case */
    function k2cCase(kebab)
    {
      return kebab.replace(/(-)([a-z])/g, (ignore, dash, char) => char.toUpperCase());
    }

    if (_dcpConfig && _dcpConfig.scheduler && _dcpConfig.scheduler.location && _dcpConfig.scheduler.location.href)
      schedulerURL = new URL(_dcpConfig.scheduler.location.href);
    else if (thisScript.getAttribute('scheduler'))
      schedulerURL = new URL(thisScript.getAttribute('scheduler'));

    if (!dcpConfigHref) {
      if (schedulerURL)
        dcpConfigHref = schedulerURL.origin + schedulerURL.pathname + 'etc/dcp-config.js' + (schedulerURL.search || thisScriptURL.search);
      else
        dcpConfigHref = thisScriptURL.origin + thisScriptURL.pathname.replace(/\/dcp-client\/dcp-client.js$/, '/etc/dcp-config.js') + thisScriptURL.search;
    }

    /**
     * <script src="dcp-client" dcp-config.scheduler.task-distributor.setting="string:foo">
     * Supported "types":
     * - json:
     * - url:
     * - string:
     * - number:
     * - true
     * - false
     */
    for (let attrib of thisScript.attributes)
    {
      let node;
      let { name, value } = attrib;
      if (!name.startsWith('dcp-config.'))
        continue;
      if (value.startsWith('string:'))
        value = value.slice(7);
      else if (value.startsWith('number:'))
        value = Number(value.slice(7));
      else if (value.startsWith('json:'))
        value = JSON.parse(value.slice(5));
      else if (value === 'true')
        value = true;
      else if (value === 'false')
        value = false;
      else if (value.startsWith('url:'))
        value = new URL(value.slice(4));
      else
      {
        console.error(`dcp-client: invalid value for attribute ${name}`);
        continue;
      }

      let edge;
      let path = name.split('.');
      path.shift();
      if (!_dcpConfig)
        _dcpConfig = {};
      node = _dcpConfig;
      while (path.length > 1)
      {
        edge = k2cCase(path.shift());
        if (!node.hasOwnProperty(edge))
          node[edge] = {};
        node = node[edge];
      }
      edge = k2cCase(path.shift());
      node[edge] = value;
    }

    /** Load dcp-config.js from scheduler, and merge with running dcpConfig */
    function loadConfig()
    {
      configScript = document.createElement('SCRIPT');
      configScript.setAttribute('type', 'text/javascript');
      configScript.setAttribute('src', dcpConfigHref);
      configScript.setAttribute('id', '_dcp_config');

      if (!thisScript.id)
        thisScript.id = '_dcp_client_loader';

      let html = configScript.outerHTML;

      /* If we know about an alternate scheduler location - by any means but usually script attribute -
       * add it into our local configuration delta; generate this delta as-needed.
       */
      if (schedulerURL)
      {
        if (!_dcpConfig)
          _dcpConfig = {};
        if (!_dcpConfig.scheduler)
          _dcpConfig.scheduler = {};
        _dcpConfig.scheduler.location = schedulerURL;
      }

      /* Preserve the config delta so that it can be merged on top of the remote config, before the
       * bundle is completely initialized. The global dcpConfig is replaced by this new script.
       */
      if (_dcpConfig)
      {
        thisScript.mergeConfig = mergeConfig;
        html += `<script>document.getElementById('${thisScript.id}').mergeConfig();</scr` + 'ipt>';
      }

      document.write(html);
      configScript.onerror = (function(e) {
        alert('Error DCP-1001: Could not load or parse scheduler configuration from URL ("' + configScript.getAttribute('src') + '")');
        console.error('dcpConfig load error: ', e);
      }).toString();
    }

    function mergeConfig()
    {
      const mergedConf = leafMerge(dcpConfig, _dcpConfig);
      Object.assign(dcpConfig, mergedConf);

      function leafMerge() /* lifted from dcp-client obj-merge.js c32e780fae88071df1bb4aebe3282220d518260e */
      {
        var target = {};

        for (let i=0; i < arguments.length; i++)
        {
          let neo = arguments[i];
          if (neo === undefined)
            continue;

          for (let prop in neo)
          {
            if (!neo.hasOwnProperty(prop))
              continue;
            if (typeof neo[prop] === 'object' && neo[prop] !== null && !Array.isArray(neo[prop]) && ['Function','Object'].includes(neo[prop].constructor.name))
              target[prop] = leafMerge(target[prop], neo[prop]);
            else
              target[prop] = neo[prop];
          }
        }

        return target;
      }
    }

    /**
     * Utility code for loading content relative to the Distributive CDN or Auth subsystems, using
     * information about their locations gleaned from dcpConfig. Since the browser is unable to preload
     * content based on these tags, page load will be slightly slower than using hard-coded tags. Note
     * also that loading scripts this way has may cause them to load after the DOMContentLoadeded event
     * has fired; caveat elitor.
     *
     * Examples:
     * 1. load stylesheet relative to CDN
     *    <link type="text/css" type="stylesheet" dcp-cdn-href="/path/to/dcp-style.css">
     * 2. load javascript relative to auth host
     *    <script dcp-auth-src="/path/to/dcp-script.js"></script>
     */
    function initMagicAttribs()
    {
      if (true
          && thisScript.getAttribute('disable-magic-attribs')
          && thisScript.getAttribute('disable-magic-attribs') !== 'false')
        return;

      function resolveWithQueryString(dcpUrl, urlTail)
      {
        const components = urlTail.split('?');
        const pathname = components[0];
        const queryString = components[1] ? ('?' + components[1]) : '';
        return String(dcpUrl.resolve(pathname)) + queryString;
      }

      function scanForTags()
      {
        if (scanForTags.busy)
          return;
        scanForTags.busy = true;
        try
        {
          const magic = { cdn: 'cdn', auth: 'pxAuth' }; /* xlate attrib fragment to conf prop */
          for (const service in magic)
          {
            const srcAttribName  = `dcp-${service}-src`;
            const hrefAttribName = `dcp-${service}-href`;
            const scripts = document.querySelectorAll(`SCRIPT[${srcAttribName}]`);
            const location = dcpConfig[magic[service]].location;

            for (const script of scripts)
            {
              const src = script.getAttribute(srcAttribName);
              script.removeAttribute(srcAttribName);
              script.src = resolveWithQueryString(location, src);
            }

            const links = document.querySelectorAll(`LINK[${hrefAttribName}]`);
            for (const link of links)
            {
              const href = link.getAttribute(hrefAttribName);
              link.removeAttribute(hrefAttribName);
              link.href = resolveWithQueryString(location, href);
            }
          }
        }
        catch (error)
        {
          console.error(error);
          throw error;
        }
        finally
        {
          scanForTags.busy = false;
        }
      } /* scanForTags */

      document.addEventListener('readystatechange', scanForTags);

      /* Watch DOM for mutations, initialize new magic relative attrib tags as they appear. */
      function setupMutationObserver()
      {
        function handleMutation(mutationList, observer)
        {
          for (const mutation of mutationList)
          {
            if (true
                && mutation.type === 'childList'
                && mutation.target.querySelectorAll('LINK:not([href]), SCRIPT:not([src])').length)
            {
              scanForTags();
              break;
            }
          }
        }

        scanForTags(); /* handle mutations that were missed to due startup races */
        const observer = new MutationObserver(handleMutation);
        observer.observe(document.head, { childList: true, subtree: true });

        if (document.body)
          observer.observe(document.body, { childList: true, subtree: true });
        else
          window.addEventListener('DOMContentLoaded', () => observer.observe(document.body, { childList: true, subtree: true }));
      }
    }

    /**
     * This function is never run directly; it is stringified and emitted in a
     * SCRIPT tag that is injected into the document. As such, it cannot close
     * over any non-global variables.
     */
    async function bundleReadyIIFE() {
      const bundleScript = document.getElementById("_dcp_client_bundle");
      const ready        = bundleScript.getAttribute('onready');
      const dcp          = bundleScript.exports;
      const KVIN         = new dcp.kvin.KVIN();
      const thisScript   = document.currentScript;
      const ihc          = (true
                            && thisScript.hasAttribute('inherit-config')
                            && thisScript.getAttribute('inherit-config') !== 'false'
                            && (window !== top || window.opener));

      KVIN.userCtors.dcpUrl$$DcpURL  = dcp['dcp-url'].DcpURL;
      KVIN.userCtors.dcpEth$$Address = dcp.wallet.Address;

      if (typeof module !== 'undefined' && typeof module.declare !== 'undefined')
        require('/internal/dcp/cjs2-shim').init(bundleScript.exports); /* CommonJS Modules/2.0d8 environment (BravoJS, NobleJS) */
      else
        window.dcp = dcp; /* vanilla JS */

      /** Let protocol know where we got out config from, so origin can be reasoned about vis a vis security */
      if (ihc)
        dcp.protocol.setSchedulerConfigLocation_fromScript((window.opener || top).document.getElementById("_dcp_config"));
      else
        dcp.protocol.setSchedulerConfigLocation_fromScript(document.getElementById("_dcp_config"));

      /**
       * Slide baked-in config underneath the remote config to provide default values, unless inheriting
       * config, in which case we deep-clone to get correct constructors for this context etc
       */
      if (!ihc)
      {
        dcpConfig = dcp['dcp-config'] = dcp.utils.leafMerge(KVIN.unmarshal(dcp['dcp-default-config']), dcp['dcp-config']);
        /**
         * Transform instances of Address-like values into Addresses. Necessary since
         * the config can't access the Address class before the bundle is loaded.
         */
        dcp.wallet.Address.patchUp(dcpConfig);
        dcp['dcp-url'].patchup(dcpConfig);
      }
      else
      {
        const parent = window.opener || top;
        dcpConfig = dcp['dcp-config'] = KVIN.unmarshal(await parent.dcp.kvin.marshalAsync(parent.dcpConfig));
      }

      if (ready)
        window.setTimeout(function bundleReadyFire() { let indirectEval=eval; indirectEval(ready) }, 0);
      window.dispatchEvent(new CustomEvent('dcpclientbundleready', { detail: { dcp, KVIN } })); /* fires as soon as bundle has been parsed */
      window.dispatchEvent(new CustomEvent('dcpclientready', { detail: dcp })); /* fires as soon as bundle has been parsed and dcp-config is ready */
    }

    /* Load dcp-client bundle from the same location as this module, extract the exports
     * from it, and attach them to the global dcp object.
     *
     * @param {boolean}   ihc           true if we should inject the inherit-config attribute onto the
     *                                  bundle's script element
     */
    function loadBundle(ihc) {
      var bundleScript = document.createElement('SCRIPT');
      var bundleSrc = thisScript.getAttribute("bundle") || (thisScript.src.replace('/dcp-client.js', '/dist/dcp-client-bundle.js'));

      bundleScript.setAttribute('type', 'text/javascript');
      bundleScript.setAttribute('src', bundleSrc);
      bundleScript.setAttribute('id', '_dcp_client_bundle');
      bundleScript.setAttribute('dcp-env', 'vanilla-web');
      bundleScript.setAttribute('onerror', `alert('Error DCP-1002: Could not load dcp-client bundle from URL ("${bundleSrc}")')`);
      bundleScript.setAttribute('onload', thisScript.getAttribute('onload'));
      thisScript.removeAttribute('onload');
      bundleScript.setAttribute('onready', thisScript.getAttribute('onready'));
      thisScript.removeAttribute('onready');
      document.write(bundleScript.outerHTML);
      document.write(`<script id='_dcp_bundleReadyIIFE' inherit-config="${ihc}">/* bundleReadyIIFE */;(${bundleReadyIIFE})()</scr` + `ipt>`);
      bundleScript = document.getElementById('_dcp_client_bundle');
      if (bundleScript)
        bundleScript.onerror = function(e) {
          console.error('Bundle load error:', e);
          bundleScript.removeAttribute('onready');
        };
    }

    /* Load the default CSS. This is enough to make the modals, etc, work without affecting the user
     * program appearance. If desired, a full DCP style can be loaded from the Distributive CDN; i.e.
     * <link rel="stylesheet" type="text/css" dcp-cdn-href="/css/dcp-style.css">
     */
    function loadCSS () {
      var href = thisScript.getAttribute('load-css');
      if (href === 'false')
        return;
      href = thisScriptURL.origin + thisScriptURL.pathname
        .replace(/\/dcp-client\/dcp-client.js$/, '/dcp-client/assets/dcp-client.css');
      const head = document.getElementsByTagName('head')[0];
      let styleLink = document.createElement('link');
      styleLink.rel = 'stylesheet';
      styleLink.href = href;
      head.appendChild(styleLink);
    }

    /**
     * Determine if an inherit attribute has been set on an element IFF we are in a heritable context;
     * eg the document has been loaded into a popup or frame with the same origin.
     */
    function hasInheritAttribute(element, suffix)
    {
      return (true
              && element.hasAttribute(`inherit-${suffix}`)
              && element.getAttribute(`inherit-${suffix}`) !== 'false'
              && (window !== top || window.opener)
              && window.location.origin == (window.opener || top).location.origin);
    }

    const ihc = hasInheritAttribute(thisScript, 'config');
    if (!ihc)
      loadConfig();
    else
      window.dcpConfig = (window.opener || top).dcpConfig; /* temp plumbing, re-written in bundleReadyIIFE */

    initMagicAttribs();
    loadCSS();

    if (hasInheritAttribute(thisScript, 'identity'))
    {
      window.addEventListener('dcpclientbundleready', ev => {
        const { dcp, KVIN } = ev.detail;
        const inheritedIdentity = (window.opener || top).dcp.identity.get();
        dcp.identity.set(inheritedIdentity);
      });
    }

    /* warning: bundle inheritance has context problems; avoid using this feature for non-trivial tasks /wg sep 2025 */
    if (hasInheritAttribute(thisScript, 'bundle'))
      window.dcp = window.opener ? window.opener.dcp : top.dcp;
    else
      loadBundle(ihc);
  }
}) /* namespaceIIFE */();
